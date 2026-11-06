//! Isostatic unilateral supports of a complete physical nullspace.
//!
//! With K N = 0, K q = f - B R requires D R = N^T f, D = N^T B.
//! A square nonsingular D fixes the reactions without contact iteration.
//! The original gauged Schur owner supplies deformation; D^T a then places
//! every contact at its original constitutive closure through q <- q + N a.
//! No stiffness is added, and the original force equations decide publication.
use super::supported::{closure_for_reaction, projection};
use super::*;
use fs_la::factor::{Lu, lu};

pub(in super::super) fn initialize(
    network: &mut CoupledModalSystem,
    external: &[f64],
    contacts: &[(ModalContact, ModalContactConfig)],
    nullspace: &[Vec<f64>],
    config: MultiContactConfig,
    gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    poll(Some(gate))?;
    network.require_static_initialization(external)?;
    let p = contacts.len();
    if p != nullspace.len() {
        return Err(invalid(
            "isostatic support needs exactly one contact per complete null direction; redundant support requires a load-sharing solve",
        ));
    }
    admit_setup(network, p, config)?;
    if p == 1 {
        return supported::initialize(
            network,
            external,
            &contacts[0].0,
            contacts[0].1,
            &nullspace[0],
            gate,
        );
    }
    // Additional small support-map factor, condition estimate, both solves,
    // and modal contractions. The original Schur owner admits its own work.
    let (n, count) = (network.mode_count() as u128, p as u128);
    let terms = n * count * (count + 8) + 2 * count * count * count + 32 * count * count;
    if terms > config.max_setup_terms as u128 {
        return Err(invalid(
            "isostatic support preparation exceeds max_setup_terms",
        ));
    }
    let response = StaticResponse::new_with_gauges(network, nullspace, gate)?;
    let mut points = Vec::with_capacity(p);
    for (contact, limits) in contacts {
        poll(Some(gate))?;
        let column = contact_column(network, contact, *limits)?;
        let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
            .map_err(ModalCouplingError::ContactLaw)?;
        points.push(Point {
            column,
            storage,
            config: *limits,
            gap: contact.law.gaps()[0],
        });
    }
    let mut matrix = vec![0.0; p * p];
    let mut loads = Vec::with_capacity(p);
    for (i, direction) in nullspace.iter().enumerate() {
        poll(Some(gate))?;
        loads.push(projection(direction, external)?.0);
        for (j, point) in points.iter().enumerate() {
            matrix[i * p + j] = projection(direction, &point.column)?.0;
        }
    }
    let map = SupportMap::new(matrix, p, gate)?;
    let reactions = map.solve(&loads, false, network.config.solve_relative_tolerance, gate)?;
    for (&reaction, point) in reactions.iter().zip(&points) {
        // Never turn a tensile reaction into a zero-force contact. At a
        // numerically unresolved sign boundary the conservative answer is refusal.
        if reaction < 0.0 {
            return Err(invalid(
                "an isostatic unilateral support requires a tensile reaction",
            ));
        }
        limit(
            "static contact force",
            reaction,
            point.config.maximum_force_n,
        )?;
    }

    let roundoff_factor = network
        .config
        .solve_relative_tolerance
        .min(64.0 * f64::EPSILON * (network.mode_count() + p) as f64);
    let mut forces = external.to_vec();
    let mut scale: Vec<f64> = external.iter().map(|f| f.abs()).collect();
    for (point, &reaction) in points.iter().zip(&reactions) {
        poll(Some(gate))?;
        for j in 0..forces.len() {
            let mapped = finite(point.column[j] * reaction)?;
            if point.column[j] != 0.0 && reaction != 0.0 && mapped == 0.0 {
                return Err(invalid(
                    "isostatic support force is below representable range",
                ));
            }
            // Keep every original subtraction; the allowance covers only
            // arithmetic cancellation of the original force terms.
            forces[j] = finite(forces[j] - mapped)?;
            scale[j] = finite(scale[j] + mapped.abs())?;
        }
    }
    let roundoff: Vec<f64> = scale
        .iter()
        .map(|s| finite(roundoff_factor * s))
        .collect::<Result<_, _>>()?;
    let mut q =
        response.solve_with_roundoff(network, &forces, true, true, Some(&roundoff), gate)?;
    let mut targets = Vec::with_capacity(p);
    for (point, &reaction) in points.iter().zip(&reactions) {
        poll(Some(gate))?;
        targets.push(if reaction == 0.0 {
            // D^T is invertible, so these independently chosen nonpenetrating
            // closures are jointly attainable. Preserve each open closure;
            // remove only interference, rather than prescribing a fake preload.
            dot(&point.column, &q)?.min(point.gap)
        } else {
            closure_for_reaction(point, reaction, gate)?.0
        });
    }
    let mut sweeps = 0;
    for sweep in 1..=config.max_sweeps {
        poll(Some(gate))?;
        sweeps = sweep;
        let differences: Vec<f64> = points
            .iter()
            .zip(&targets)
            .map(|(point, &target)| finite(target - dot(&point.column, &q)?))
            .collect::<Result<_, _>>()?;
        let translation = map.solve(
            &differences,
            true,
            network.config.solve_relative_tolerance,
            gate,
        )?;
        for (direction, &distance) in nullspace.iter().zip(&translation) {
            poll(Some(gate))?;
            for (q, &entry) in q.iter_mut().zip(direction) {
                *q = finite(*q + entry * distance)?;
            }
        }
        let mut resolved = true;
        for (point, &reaction) in points.iter().zip(&reactions) {
            poll(Some(gate))?;
            let expected = point.force(dot(&point.column, &q)?)?;
            resolved &= finite(reaction - expected)?.abs()
                <= force_tolerance(reaction, expected, point.config)?;
        }
        if resolved || translation.iter().all(|a| *a == 0.0) {
            break;
        }
    }
    // The small-map equations never replace actual original spring, modal,
    // contact, material-history and stored-energy checks on the restored state.
    let energy =
        network.stage_static_equilibrium_with_allowance(&q, &forces, Some(&roundoff), gate)?;
    publish(
        network, external, &points, None, &reactions, energy, sweeps, gate,
    )
}

// The original fs-la pivoted factor serves force balance and its transpose
// closure solve. Row/column equilibration prevents the physical mass/shape
// units from deciding rank. The condition screen is a floating-point estimate,
// not an interval rank or mechanical stability certificate.
struct SupportMap {
    matrix: Vec<f64>,
    scaled: Vec<f64>,
    row: Vec<f64>,
    column: Vec<f64>,
    factor: Lu,
}
impl SupportMap {
    fn new(matrix: Vec<f64>, p: usize, gate: &CancelGate) -> Result<Self, ModalCouplingError> {
        let row: Vec<f64> = matrix
            .chunks_exact(p)
            .map(|r| r.iter().map(|v| v.abs()).fold(0.0, f64::max))
            .collect();
        if row.contains(&0.0) {
            return Err(invalid(
                "isostatic contacts do not restrain every declared free direction",
            ));
        }
        let mut scaled = Vec::with_capacity(matrix.len());
        for (i, entries) in matrix.chunks_exact(p).enumerate() {
            poll(Some(gate))?;
            for &entry in entries {
                scaled.push(divide(entry, row[i])?);
            }
        }
        let mut column = vec![0.0_f64; p];
        for entries in scaled.chunks_exact(p) {
            for (j, &entry) in entries.iter().enumerate() {
                column[j] = column[j].max(entry.abs());
            }
        }
        if column.contains(&0.0) {
            return Err(invalid(
                "an isostatic contact is blind to the complete free subspace",
            ));
        }
        for entries in scaled.chunks_exact_mut(p) {
            for (entry, &scale) in entries.iter_mut().zip(&column) {
                *entry = divide(*entry, scale)?;
            }
        }
        poll(Some(gate))?;
        let factor = lu(&scaled, p).map_err(ModalCouplingError::Factor)?;
        let condition = finite(factor.condition_1(&scaled))?;
        if condition * 256.0 * f64::EPSILON * p as f64 >= 1.0 {
            return Err(invalid(
                "isostatic contact directions are dependent or unresolved at roundoff",
            ));
        }
        poll(Some(gate))?;
        Ok(Self {
            matrix,
            scaled,
            row,
            column,
            factor,
        })
    }

    fn solve(
        &self,
        rhs: &[f64],
        transpose: bool,
        tolerance: f64,
        gate: &CancelGate,
    ) -> Result<Vec<f64>, ModalCouplingError> {
        poll(Some(gate))?;
        let (before, after) = if transpose {
            (&self.column, &self.row)
        } else {
            (&self.row, &self.column)
        };
        let b: Vec<f64> = rhs
            .iter()
            .zip(before)
            .map(|(&r, &s)| divide(r, s))
            .collect::<Result<_, _>>()?;
        let mut x = b.clone();
        if transpose {
            self.factor.solve_transpose(&mut x);
        } else {
            self.factor.solve(&mut x);
        }
        residual(&self.scaled, &x, &b, transpose, tolerance)?;
        for (value, &scale) in x.iter_mut().zip(after) {
            *value = divide(*value, scale)?;
        }
        residual(&self.matrix, &x, rhs, transpose, tolerance)?;
        poll(Some(gate))?;
        Ok(x)
    }
}

fn divide(value: f64, scale: f64) -> Result<f64, ModalCouplingError> {
    let result = finite(value / scale)?;
    if value != 0.0 && result == 0.0 {
        return Err(invalid(
            "isostatic support scaling is below representable range",
        ));
    }
    Ok(result)
}

fn residual(
    matrix: &[f64],
    x: &[f64],
    rhs: &[f64],
    transpose: bool,
    tolerance: f64,
) -> Result<(), ModalCouplingError> {
    let p = rhs.len();
    let mut relative = 0.0_f64;
    for i in 0..p {
        let mut applied = 0.0;
        let mut scale = rhs[i].abs();
        for j in 0..p {
            let a = matrix[if transpose { j * p + i } else { i * p + j }];
            let term = finite(a * x[j])?;
            if a != 0.0 && x[j] != 0.0 && term == 0.0 {
                return Err(invalid(
                    "isostatic support residual is below representable range",
                ));
            }
            applied = finite(applied + term)?;
            scale = finite(scale + term.abs())?;
        }
        let defect = finite(rhs[i] - applied)?.abs();
        if scale > 0.0 {
            relative = relative.max(defect / scale);
        }
    }
    if relative > tolerance {
        return Err(ModalCouplingError::SolveResidual {
            relative,
            tolerance,
        });
    }
    Ok(())
}
