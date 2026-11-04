//! Static reciprocal maps between original physical power-port banks.
//!
//! This is a composition of the original state, storage and losses. The native
//! integrator evolves its complete coupled operator; there is no held-interface
//! approximation, extra state, rank cutoff, or independent coupling integrator.

use crate::{PhsError, PortHamiltonian, SumStorage, compose_resistance_witness};

/// Explicit capacity for a complete dense physical interconnection.
///
/// These bounds cover composition, not subsequent time integration. Original
/// source owners are already allocated. Native skew/PSD and resistance-witness
/// admission still run, with their work bounded by `max_states`.
#[derive(Clone, Copy, Debug)]
pub struct InterconnectionBudget {
    /// Maximum combined original state dimension, including every memory/mode.
    pub max_states: usize,
    /// Maximum combined ORIGINAL port count, including ports closed internally.
    pub max_total_ports: usize,
    /// Conservative additional dense storage: 64*n*n + n*external_ports +
    /// twice the contraction scratch extent in doubles, plus one byte per original port. Includes
    /// native matrix/factor admission scratch, but not source-owned storage.
    pub max_dense_bytes: usize,
    /// Maximum product/FMA terms in the two coupling contractions. Each product
    /// retains its FMA residual, and both first-stage words enter stage two.
    /// Zero map entries contribute no terms; nonzero coefficients are never cut.
    pub max_coupling_multiply_terms: u64,
}

/// Compose the ORIGINAL systems through a reciprocal rectangular physical map.
///
/// With outputs in the selected bank orders, `u_right = M*y_left` and
/// `u_left = -M^T*y_right`. `right_from_left` is row-major RIGHT-by-LEFT.
/// Selected ports become internal; unselected left ports remain external in
/// original order, followed by unselected right ports in original order.
/// States remain `[left_state, right_state]`, including every null direction.
/// An all-zero map closes the selected inputs at zero without removing states.
/// Empty banks (both empty) give the ordinary direct sum.
///
/// Piston areas, orientations, transformer ratios and already assembled spatial
/// weights are physical coefficients. Port coordinates must already be dual
/// with respect to power, including quadrature/integration measures. This API
/// does not construct or certify a spatial projection. Coefficients are fixed
/// in the returned system; work from moving geometry or changed storage needs
/// its own original owner. There are no numerical input scales or port cutoff.
///
/// The full block-skew operator is passed to the existing pHS owner. Quadratic
/// source storage therefore supports the existing exact ZOH/tone preparation;
/// nonlinear source storage retains its original discrete-gradient law.
/// No new exponential, integrator, damping, or storage law is introduced here.
///
/// # Errors
/// Refuses duplicate/out-of-range ports, malformed/nonfinite data, exceeded
/// explicit capacities, unrepresentable contraction arithmetic, or original
/// structural/resistance admission. Inputs are consumed only into this result.
pub fn interconnect_mapped(
    left: PortHamiltonian,
    right: PortHamiltonian,
    left_ports: &[usize],
    right_ports: &[usize],
    right_from_left: &[f64],
    budget: InterconnectionBudget,
) -> Result<PortHamiltonian, PhsError> {
    assemble(
        left,
        right,
        Connection::Mapped {
            left: left_ports,
            right: right_ports,
            coefficients: right_from_left,
        },
        Some(budget),
    )
}

pub(super) fn paired(
    left: PortHamiltonian,
    right: PortHamiltonian,
    pairs: &[(usize, usize)],
) -> Result<PortHamiltonian, PhsError> {
    assemble(left, right, Connection::Paired(pairs), None)
}

pub(super) fn transformer(
    left: PortHamiltonian,
    right: PortHamiltonian,
    port_left: usize,
    port_right: usize,
    ratio: f64,
) -> Result<PortHamiltonian, PhsError> {
    // The legacy transformer has the opposite sign to the canonical map.
    if port_left >= left.m || port_right >= right.m {
        return Err(PhsError::BadPortPairing);
    }
    if !ratio.is_finite() {
        return Err(refuse("transformer ratio"));
    }
    assemble(
        left,
        right,
        Connection::Transformer {
            left: port_left,
            right: port_right,
            ratio,
        },
        None,
    )
}

enum Connection<'a> {
    Paired(&'a [(usize, usize)]),
    Transformer {
        left: usize,
        right: usize,
        ratio: f64,
    },
    Mapped {
        left: &'a [usize],
        right: &'a [usize],
        coefficients: &'a [f64],
    },
}

fn bad(what: &'static str) -> PhsError {
    PhsError::Dimension { what }
}
fn refuse(what: &'static str) -> PhsError {
    PhsError::Unadmitted { what }
}
fn filled<T: Clone>(count: usize, value: T) -> Result<Vec<T>, PhsError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| bad("interconnection allocation capacity"))?;
    result.resize(count, value);
    Ok(result)
}
fn source_extent(system: &PortHamiltonian) -> Result<(), PhsError> {
    if Some(system.j.len()) != system.n.checked_mul(system.n)
        || system.r.len() != system.j.len()
        || Some(system.g.len()) != system.n.checked_mul(system.m)
    {
        return Err(bad("interconnection original matrix extents"));
    }
    if system
        .j
        .iter()
        .chain(&system.r)
        .chain(&system.g)
        .any(|value| !value.is_finite())
    {
        return Err(refuse("interconnection requires finite original matrices"));
    }
    Ok(())
}
fn mark(mask: &mut [bool], port: usize) -> Result<(), PhsError> {
    if port >= mask.len() || mask[port] {
        return Err(PhsError::BadPortPairing);
    }
    mask[port] = true;
    Ok(())
}

#[derive(Clone, Copy)]
struct Contraction {
    left_first: bool,
    scratch: usize,
    terms: u128,
}
fn contraction(
    na: usize,
    nb: usize,
    nl: usize,
    nr: usize,
    nonzero: usize,
) -> Result<Contraction, PhsError> {
    if nonzero == 0 || na == 0 || nb == 0 {
        return Ok(Contraction {
            left_first: true,
            scratch: 0,
            terms: 0,
        });
    }
    let cost = |rows: usize, other: usize, columns: usize| {
        (rows as u128)
            .checked_mul(nonzero as u128)
            .and_then(|terms| terms.checked_mul(2))
            .and_then(|first| {
                (rows as u128)
                    .checked_mul(other as u128)
                    .and_then(|cross| cross.checked_mul(columns as u128))
                    .and_then(|terms| terms.checked_mul(4))
                    .and_then(|second| first.checked_add(second))
            })
            .ok_or_else(|| bad("interconnection contraction work extent"))
    };
    let left_cost = cost(na, nb, nr)?;
    let right_cost = cost(nb, na, nl)?;
    let left_first = left_cost <= right_cost;
    let scratch = if left_first {
        na.checked_mul(nr)
    } else {
        nb.checked_mul(nl)
    }
    .ok_or_else(|| bad("interconnection contraction scratch extent"))?;
    Ok(Contraction {
        left_first,
        scratch,
        terms: left_cost.min(right_cost),
    })
}

mod accumulation;
use accumulation::Sum;

fn checked<T>(value: Option<T>) -> Result<T, PhsError> {
    value.ok_or_else(|| refuse("unrepresentable interconnection contraction"))
}

fn mapped_cross(
    a: &PortHamiltonian,
    b: &PortHamiltonian,
    left: &[usize],
    right: &[usize],
    coefficients: &[f64],
    plan: Contraction,
    j: &mut [f64],
) -> Result<(), PhsError> {
    if plan.terms == 0 {
        return Ok(());
    }
    let (na, nb, nl, nr) = (a.n, b.n, left.len(), right.len());
    let n = na + nb;
    // Do not round a first-stage dot product into one f64: cancellation between
    // second-stage columns can make its low word the entire physical coupling.
    let mut temporary = filled(plan.scratch, Sum::default())?;
    if plan.left_first {
        // (G_left M^T) G_right^T; O(na*nnz(M) + na*nb*nr) products.
        for i in 0..na {
            for r in 0..nr {
                let mut value = Sum::default();
                for l in 0..nl {
                    let coefficient = coefficients[r * nl + l];
                    if coefficient != 0.0 {
                        checked(value.product(a.g[i * a.m + left[l]], coefficient))?;
                    }
                }
                temporary[i * nr + r] = value;
            }
        }
        for i in 0..na {
            for k in 0..nb {
                let mut value = Sum::default();
                for r in 0..nr {
                    checked(value.add_scaled(temporary[i * nr + r], b.g[k * b.m + right[r]]))?;
                }
                let value = checked(value.value())?;
                j[i * n + na + k] = -value;
                j[(na + k) * n + i] = value;
            }
        }
    } else {
        // G_left (G_right M)^T; choose this for the smaller native contraction.
        for k in 0..nb {
            for l in 0..nl {
                let mut value = Sum::default();
                for r in 0..nr {
                    let coefficient = coefficients[r * nl + l];
                    if coefficient != 0.0 {
                        checked(value.product(b.g[k * b.m + right[r]], coefficient))?;
                    }
                }
                temporary[k * nl + l] = value;
            }
        }
        for i in 0..na {
            for k in 0..nb {
                let mut value = Sum::default();
                for l in 0..nl {
                    checked(value.add_scaled(temporary[k * nl + l], a.g[i * a.m + left[l]]))?;
                }
                let value = checked(value.value())?;
                j[i * n + na + k] = -value;
                j[(na + k) * n + i] = value;
            }
        }
    }
    Ok(())
}

fn assemble(
    a: PortHamiltonian,
    b: PortHamiltonian,
    connection: Connection<'_>,
    budget: Option<InterconnectionBudget>,
) -> Result<PortHamiltonian, PhsError> {
    source_extent(&a)?;
    source_extent(&b)?;
    let (na, nb) = (a.n, b.n);
    let n = na
        .checked_add(nb)
        .ok_or_else(|| bad("interconnection state extent"))?;
    let square = n
        .checked_mul(n)
        .ok_or_else(|| bad("interconnection matrix extent"))?;
    let total_ports =
        a.m.checked_add(b.m)
            .ok_or_else(|| bad("interconnection port extent"))?;
    if let Some(budget) = budget {
        if n > budget.max_states || total_ports > budget.max_total_ports {
            return Err(refuse("interconnection state or port capacity exceeded"));
        }
    }
    let (nl, nr, plan) = match &connection {
        Connection::Paired(pairs) => (pairs.len(), pairs.len(), None),
        Connection::Transformer { .. } => (1, 1, None),
        Connection::Mapped {
            left,
            right,
            coefficients,
        } => {
            if left.is_empty() != right.is_empty()
                || Some(coefficients.len()) != left.len().checked_mul(right.len())
                || coefficients.iter().any(|value| !value.is_finite())
            {
                return Err(bad(
                    "finite complete rectangular interconnection map required",
                ));
            }
            let nonzero = coefficients.iter().filter(|&&value| value != 0.0).count();
            (
                left.len(),
                right.len(),
                Some(contraction(na, nb, left.len(), right.len(), nonzero)?),
            )
        }
    };
    if nl > a.m || nr > b.m {
        return Err(PhsError::BadPortPairing);
    }
    let m = (a.m - nl)
        .checked_add(b.m - nr)
        .ok_or_else(|| bad("interconnection external port extent"))?;
    let input_entries = n
        .checked_mul(m)
        .ok_or_else(|| bad("interconnection external source extent"))?;
    if let Some(budget) = budget {
        let plan = plan.ok_or_else(|| bad("mapped interconnection contraction plan"))?;
        let bytes = (64 * square as u128 + input_entries as u128 + 2 * plan.scratch as u128)
            * core::mem::size_of::<f64>() as u128
            + total_ports as u128;
        if bytes > budget.max_dense_bytes as u128 {
            return Err(refuse("interconnection dense storage capacity exceeded"));
        }
        if plan.terms > budget.max_coupling_multiply_terms as u128 {
            return Err(refuse("interconnection coupling work capacity exceeded"));
        }
    }
    let mut used_a = filled(a.m, false)?;
    let mut used_b = filled(b.m, false)?;
    match &connection {
        Connection::Paired(pairs) => {
            for &(left, right) in *pairs {
                mark(&mut used_a, left)?;
                mark(&mut used_b, right)?;
            }
        }
        Connection::Transformer { left, right, .. } => {
            mark(&mut used_a, *left)?;
            mark(&mut used_b, *right)?;
        }
        Connection::Mapped { left, right, .. } => {
            for &port in *left {
                mark(&mut used_a, port)?;
            }
            for &port in *right {
                mark(&mut used_b, port)?;
            }
        }
    }
    let mut j = filled(square, 0.0)?;
    let mut r = filled(square, 0.0)?;
    let mut g = filled(input_entries, 0.0)?;
    for (system, offset) in [(&a, 0), (&b, na)] {
        for i in 0..system.n {
            let target = (offset + i) * n + offset;
            j[target..target + system.n]
                .copy_from_slice(&system.j[i * system.n..(i + 1) * system.n]);
            r[target..target + system.n]
                .copy_from_slice(&system.r[i * system.n..(i + 1) * system.n]);
        }
    }
    match connection {
        Connection::Paired(pairs) => {
            // Keep the legacy sparse pairing cost and arithmetic order.
            for &(left, right) in pairs {
                for i in 0..na {
                    for k in 0..nb {
                        let value = a.g[i * a.m + left] * b.g[k * b.m + right];
                        j[i * n + na + k] -= value;
                        j[(na + k) * n + i] += value;
                    }
                }
            }
        }
        Connection::Transformer { left, right, ratio } => {
            // Preserve u_left = ratio*y_right, u_right = -ratio*y_left.
            for i in 0..na {
                for k in 0..nb {
                    let value = ratio * a.g[i * a.m + left] * b.g[k * b.m + right];
                    j[i * n + na + k] += value;
                    j[(na + k) * n + i] -= value;
                }
            }
        }
        Connection::Mapped {
            left,
            right,
            coefficients,
        } => {
            mapped_cross(
                &a,
                &b,
                left,
                right,
                coefficients,
                plan.ok_or_else(|| bad("mapped interconnection contraction plan"))?,
                &mut j,
            )?;
        }
    }
    if j.iter().any(|value| !value.is_finite()) {
        return Err(refuse("unrepresentable interconnection structure"));
    }
    let mut column = 0;
    for (system, offset, used) in [(&a, 0, &used_a), (&b, na, &used_b)] {
        for (port, &internal) in used.iter().enumerate() {
            if !internal {
                for i in 0..system.n {
                    g[(offset + i) * m + column] = system.g[i * system.m + port];
                }
                column += 1;
            }
        }
    }
    let resistance_witness = compose_resistance_witness(&[&a, &b], n, &r)?;
    let storage = Box::new(SumStorage {
        a: (a.storage, na),
        b: (b.storage, nb),
    });
    let mut system = PortHamiltonian::new(n, m, j, r, g, storage)?;
    system.resistance_witness = resistance_witness;
    Ok(system)
}
