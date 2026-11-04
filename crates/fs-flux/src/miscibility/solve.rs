//! Bounded scalar convex-split exchanges and one atomic conservative sweep.

use super::{
    BinaryMiscibilityCandidate, BinaryMiscibilityError as Error, BinaryMiscibilityNetwork,
    BinaryMiscibilitySolve, BinaryMiscibilityState, BinaryMiscibilityStep, Result, Sum, finite,
    nonnegative, poll, positive,
};

pub(super) fn step(
    network: &BinaryMiscibilityNetwork,
    initial: &BinaryMiscibilityState,
    duration_s: f64,
    controls: BinaryMiscibilitySolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<BinaryMiscibilityCandidate> {
    poll(&mut cancelled)?;
    network.admit(initial)?;
    if !positive(duration_s)
        || !(1..=256).contains(&controls.maximum_root_iterations)
        || !nonnegative(controls.absolute_transfer_tolerance_mol)
        || !controls.relative_transfer_tolerance.is_finite()
        || controls.relative_transfer_tolerance < 4.0 * f64::EPSILON
        || controls.relative_transfer_tolerance > 1e-3
        || !nonnegative(controls.amount_balance_tolerance_mol)
        || !nonnegative(controls.free_energy_tolerance_j)
    {
        return Err(Error::Invalid(
            "positive time and bounded root/conservation/energy controls required",
        ));
    }
    let elapsed_s = finite(initial.elapsed_s + duration_s, "elapsed time")?;
    if elapsed_s <= initial.elapsed_s {
        return Err(Error::Arithmetic("time increment is unrepresentable"));
    }
    let before = network.free_energy_j(initial, &mut cancelled)?;
    let mut state = initial.clone();
    let mut gradient = network.gradient(&state.fractions, &mut cancelled)?;
    let mut transfers = vec![0.0; network.owner.faces.len()];
    let mut iterations = 0;
    let mut maximum_residual: f64 = 0.0;
    let mut dissipation = Sum::default();
    for (face_index, face) in network.owner.faces.iter().enumerate() {
        poll(&mut cancelled)?;
        let Some(root) = solve_face(
            network,
            &state,
            &gradient,
            face_index,
            duration_s,
            controls,
            &mut cancelled,
        )?
        else {
            continue;
        };
        let [i, j] = face.physical.cells;
        if root.transfer != 0.0 {
            // Each represented debit and credit must move in the declared
            // direction. A positive incoming amount may not disappear beside
            // a large inventory while its opposite debit remains represented.
            let signs = if root.transfer > 0.0 {
                [-1.0, 1.0, 1.0, -1.0]
            } else {
                [1.0, -1.0, -1.0, 1.0]
            };
            let differences = [
                root.amounts[0][0] - state.amounts[i][0],
                root.amounts[0][1] - state.amounts[i][1],
                root.amounts[1][0] - state.amounts[j][0],
                root.amounts[1][1] - state.amounts[j][1],
            ];
            if differences
                .into_iter()
                .zip(signs)
                .any(|(delta, sign)| delta * sign <= 0.0)
            {
                return Err(Error::Arithmetic(
                    "nonzero component transfer lost in inventory arithmetic",
                ));
            }
            for (cell, fraction) in [(i, root.fractions[0]), (j, root.fractions[1])] {
                let change = fraction[0] - state.fractions[cell][0];
                gradient[cell].add(finite(
                    network.owner.gradient_diagonal_j[cell] * change,
                    "local gradient update",
                )?)?;
                let range =
                    network.owner.adjacent_offsets[cell]..network.owner.adjacent_offsets[cell + 1];
                for &adjacent in &network.owner.adjacent_faces[range] {
                    poll(&mut cancelled)?;
                    let neighbor_face = &network.owner.faces[adjacent];
                    let [a, b] = neighbor_face.physical.cells;
                    let neighbor = if a == cell { b } else { a };
                    gradient[neighbor].add(finite(
                        -neighbor_face.gradient_j * change,
                        "neighbor gradient update",
                    )?)?;
                }
            }
            state.amounts[i] = root.amounts[0];
            state.amounts[j] = root.amounts[1];
            state.amount_remainders[i] = root.amount_remainders[0];
            state.amount_remainders[j] = root.amount_remainders[1];
            state.fractions[i] = root.fractions[0];
            state.fractions[j] = root.fractions[1];
        }
        let dissipated = finite(
            root.transfer * root.chemical_drop,
            "convex-split dissipation",
        )?;
        if dissipated < 0.0 {
            return Err(Error::Arithmetic(
                "implicit exchange opposes its chemical-potential drop",
            ));
        }
        dissipation.add(dissipated)?;
        transfers[face_index] = root.transfer;
        iterations += root.iterations;
        maximum_residual = maximum_residual.max(root.residual.abs());
    }
    let mut balances = [Sum::default(); 2];
    let mut capacity_residual: f64 = 0.0;
    for (i, (&old, &new)) in initial.amounts.iter().zip(&state.amounts).enumerate() {
        poll(&mut cancelled)?;
        for species in 0..2 {
            balances[species].add(finite(new[species] - old[species], "component balance")?)?;
        }
        capacity_residual =
            capacity_residual.max((new[0] + new[1] - network.owner.capacities_mol[i]).abs());
    }
    let balances = [balances[0].value()?, balances[1].value()?];
    if balances
        .iter()
        .any(|x| x.abs() > controls.amount_balance_tolerance_mol)
        || capacity_residual > controls.amount_balance_tolerance_mol
    {
        return Err(Error::Conservation {
            component_residual_mol: balances,
            capacity_residual_mol: capacity_residual,
        });
    }
    let after = network.free_energy_j(&state, &mut cancelled)?;
    let change = network.free_energy_change(initial, &state, &mut cancelled)?;
    let dissipation = dissipation.value()?;
    let defect = finite(change + dissipation, "whole-step free-energy inequality")?;
    if change > controls.free_energy_tolerance_j || defect > controls.free_energy_tolerance_j {
        return Err(Error::FreeEnergy {
            change_j: change,
            dissipation_j: dissipation,
            tolerance_j: controls.free_energy_tolerance_j,
        });
    }
    state.elapsed_s = elapsed_s;
    poll(&mut cancelled)?;
    Ok(BinaryMiscibilityCandidate {
        state,
        step: BinaryMiscibilityStep {
            duration_s,
            transferred_a_mol: transfers,
            root_iterations: iterations,
            maximum_root_residual_mol: maximum_residual,
            component_residual_mol: balances,
            maximum_capacity_residual_mol: capacity_residual,
            free_energy_before_j: before,
            free_energy_after_j: after,
            free_energy_change_j: change,
            convex_split_dissipation_j: dissipation,
        },
    })
}

struct Root {
    transfer: f64,
    amounts: [[f64; 2]; 2],
    amount_remainders: [[f64; 2]; 2],
    fractions: [[f64; 2]; 2],
    chemical_drop: f64,
    residual: f64,
    derivative: Option<f64>,
    iterations: usize,
}

fn solve_face(
    network: &BinaryMiscibilityNetwork,
    state: &BinaryMiscibilityState,
    gradient: &[Sum],
    face_index: usize,
    duration_s: f64,
    controls: BinaryMiscibilitySolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Option<Root>> {
    let face = &network.owner.faces[face_index];
    let [i, j] = face.physical.cells;
    if face.mobility_mol2_j_s == 0.0 {
        return Ok(None);
    }
    let old = [state.fractions[i], state.fractions[j]];
    // Arithmetic face composition is the binary Maxwell-Stefan molar-frame
    // mobility: c D x_A x_B/(RT), never a signed effective Fick diffusivity.
    let mean = [
        0.5 * old[0][0] + 0.5 * old[1][0],
        0.5 * old[0][1] + 0.5 * old[1][1],
    ];
    if (old[0][0] == 0.0 && old[1][0] == 0.0) || (old[0][1] == 0.0 && old[1][1] == 0.0) {
        return Ok(None);
    }
    let mobility = face.mobility_mol2_j_s * (mean[0] * mean[1]);
    let time_mobility = mobility * duration_s;
    if !positive(mobility) || !positive(time_mobility) {
        return Err(Error::Arithmetic(
            "positive active-face mobility or time scale is unrepresentable",
        ));
    }
    let amounts = [state.amounts[i], state.amounts[j]];
    let amount_remainders = [state.amount_remainders[i], state.amount_remainders[j]];
    let capacity = [
        network.owner.capacities_mol[i],
        network.owner.capacities_mol[j],
    ];
    let original_gradient = [gradient[i].value()?, gradient[j].value()?];
    let diagonal = [
        network.owner.gradient_diagonal_j[i],
        network.owner.gradient_diagonal_j[j],
    ];
    let mut lower = -amounts[0][1].min(amounts[1][0]);
    let mut upper = amounts[0][0].min(amounts[1][1]);
    if lower >= upper {
        return Err(Error::Arithmetic(
            "active face has no representable interior exchange interval",
        ));
    }
    let sample = |transfer: f64| -> Result<Root> {
        let mut next = [[0.0; 2]; 2];
        let mut remainders = [[0.0; 2]; 2];
        for cell in 0..2 {
            for species in 0..2 {
                let delta = if cell == species { -transfer } else { transfer };
                (next[cell][species], remainders[cell][species]) = add_amount(
                    amounts[cell][species],
                    amount_remainders[cell][species],
                    delta,
                )?;
            }
        }
        if next.iter().flatten().any(|&n| !positive(n)) {
            return Err(Error::Arithmetic(
                "implicit exchange must retain an interior positive pair",
            ));
        }
        let fractions = [
            network.fractions(next[0], capacity[0])?,
            network.fractions(next[1], capacity[1])?,
        ];
        let changes = [fractions[0][0] - old[0][0], fractions[1][0] - old[1][0]];
        let mut gradient_left = Sum::default();
        gradient_left.add(original_gradient[0])?;
        gradient_left.add(finite(diagonal[0] * changes[0], "implicit left gradient")?)?;
        gradient_left.add(finite(
            -face.gradient_j * changes[1],
            "implicit right-neighbor gradient",
        )?)?;
        let mut gradient_right = Sum::default();
        gradient_right.add(original_gradient[1])?;
        gradient_right.add(finite(diagonal[1] * changes[1], "implicit right gradient")?)?;
        gradient_right.add(finite(
            -face.gradient_j * changes[0],
            "implicit left-neighbor gradient",
        )?)?;
        let mut drop = Sum::default();
        drop.add(
            network
                .model()
                .convex_split_exchange_potential_j_per_mol(fractions[0], old[0])?,
        )?;
        drop.add(
            -network
                .model()
                .convex_split_exchange_potential_j_per_mol(fractions[1], old[1])?,
        )?;
        drop.add(finite(
            gradient_left.value()? / capacity[0],
            "left molar gradient",
        )?)?;
        drop.add(finite(
            -gradient_right.value()? / capacity[1],
            "right molar gradient",
        )?)?;
        let chemical_drop = drop.value()?;
        let residual = finite(
            transfer - time_mobility * chemical_drop,
            "implicit exchange residual",
        )?;
        // Newton is only an accelerator. Overflow of a positive curvature near
        // a pure endpoint falls back to the same bounded monotone bisection.
        let derivative = network
            .model()
            .convex_curvature_j_per_mol(fractions[0])
            .ok()
            .zip(
                network
                    .model()
                    .convex_curvature_j_per_mol(fractions[1])
                    .ok(),
            )
            .and_then(|(left, right)| {
                let curvature = left / capacity[0]
                    + right / capacity[1]
                    + diagonal[0] / capacity[0] / capacity[0]
                    + diagonal[1] / capacity[1] / capacity[1]
                    + 2.0 * (face.gradient_j / capacity[0] / capacity[1]);
                let derivative = 1.0 + time_mobility * curvature;
                positive(derivative).then_some(derivative)
            });
        Ok(Root {
            transfer,
            amounts: next,
            amount_remainders: remainders,
            fractions,
            chemical_drop,
            residual,
            derivative,
            iterations: 0,
        })
    };
    // Endpoint residual signs follow the singular log limits. Evaluate only
    // interior trial amounts; a pure initial contact needs no fictitious trace.
    let mut transfer = if old.iter().flatten().all(|&x| x > 0.0) {
        0.0
    } else {
        0.5 * lower + 0.5 * upper
    };
    let mut last_residual = f64::INFINITY;
    for iteration in 1..=controls.maximum_root_iterations {
        poll(cancelled)?;
        if transfer <= lower || transfer >= upper {
            break;
        }
        let mut root = sample(transfer)?;
        root.iterations = iteration;
        last_residual = root.residual;
        let residual_tolerance = finite(
            controls.absolute_transfer_tolerance_mol
                + controls.relative_transfer_tolerance
                    * transfer
                        .abs()
                        .max((time_mobility * root.chemical_drop).abs()),
            "root residual tolerance",
        )?;
        // Convexity gives r'(q)>=1, so |r| already bounds the molar transfer
        // error. A second narrow-bracket condition would throw away Newton's
        // convergence, and would invent tiny transfers at rounded stationarity.
        if root.residual.abs() <= residual_tolerance {
            return Ok(Some(root));
        }
        if root.residual < 0.0 {
            lower = transfer;
        } else {
            upper = transfer;
        }
        let midpoint = 0.5 * lower + 0.5 * upper;
        if midpoint == lower || midpoint == upper {
            break;
        }
        let proposal = root
            .derivative
            .map(|derivative| transfer - root.residual / derivative);
        transfer = proposal
            // The residual already controls transfer error. An artificial
            // distance from the bracket endpoints would force many bisections
            // whenever a small physical timestep produces a small transfer.
            .filter(|&next| next > lower && next < upper && next != transfer)
            .unwrap_or(midpoint);
    }
    Err(Error::Root {
        face: face_index,
        iterations: controls.maximum_root_iterations,
        residual_mol: last_residual,
    })
}

/// Add a trial debit/credit to the same original two-word amount every time.
/// These low words belong to the immutable state and survive accepted sweeps;
/// they never change capacity, renormalize inventories, or relax admission.
fn add_amount(amount: f64, remainder: f64, delta: f64) -> Result<(f64, f64)> {
    let (sum, first_error) = two_sum(amount, delta)?;
    let (sum, second_error) = two_sum(sum, remainder)?;
    two_sum(
        sum,
        finite(first_error + second_error, "component amount remainder")?,
    )
}

// General TwoSum, without the magnitude ordering assumed by QuickTwoSum:
// a near-depleted inventory may be smaller than its carried rounding word.
fn two_sum(a: f64, b: f64) -> Result<(f64, f64)> {
    let sum = finite(a + b, "component amount accumulation")?;
    let represented_b = sum - a;
    let error = (a - (sum - represented_b)) + (b - represented_b);
    Ok((sum, finite(error, "component amount rounding remainder")?))
}
