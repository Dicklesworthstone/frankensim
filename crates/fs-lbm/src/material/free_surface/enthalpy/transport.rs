//! Conservative implicit donor transport over native transferred masses.
//!
//! Row i is `(m_new_i + sum(outgoing dm))*h_i - sum(incoming dm*h_donor)
//! = m_old_i*h_old_i`. The original mass balance makes this an anchored
//! M-matrix when each participating component receives old positive storage.
//! Conversion carry is an ordinary node with its original old/new mass. No
//! artificial mass, heat, concentration floor or diagonal regularizer is used.

use super::*;
use crate::conjugate::krylov::{KrylovOutcome, bicgstab_ilu0};
use crate::d3q19::freesurface3::transport::MassTransport3;
use fs_sparse::Coo;

pub(super) fn advect(
    trace: &MassTransport3,
    old_h: &BTreeMap<MassLocation3, f64>,
    controls: FreeSurfaceEnthalpyControls,
    gate: &CancelGate,
) -> Result<(BTreeMap<MassLocation3, f64>, KrylovOutcome), Error> {
    let mut keys = BTreeSet::new();
    for (&key, &mass) in trace.before.iter().chain(&trace.after) {
        if mass > 0.0 {
            keys.insert(key);
        }
    }
    for route in &trace.routes {
        keys.insert(route.from);
        keys.insert(route.to);
    }
    if keys.is_empty() || keys.len() > controls.max_cells {
        return Err(Error::Invalid("transport node budget exhausted"));
    }
    let index: BTreeMap<_, _> = keys
        .into_iter()
        .enumerate()
        .map(|(i, key)| (key, i))
        .collect();
    let n = index.len();
    let mut old = vec![0.0; n];
    let mut new = vec![0.0; n];
    let mut incoming = vec![0.0; n];
    let mut outgoing = vec![0.0; n];
    let mut reachable = vec![false; n];
    let mut downstream = vec![Vec::new(); n];
    let mut pending = Vec::new();
    let mut b = vec![0.0; n];
    let mut h = vec![0.0; n];
    for (&key, &i) in &index {
        poll(gate)?;
        old[i] = trace.before.get(&key).copied().unwrap_or(0.0);
        new[i] = trace.after.get(&key).copied().unwrap_or(0.0);
        if old[i] > 0.0 {
            h[i] = *old_h.get(&key).ok_or(Error::Invalid(
                "positive source mass needs retained enthalpy",
            ))?;
            b[i] = old[i] * h[i];
            reachable[i] = true;
            pending.push(i);
        }
    }
    for route in &trace.routes {
        let from = index[&route.from];
        let to = index[&route.to];
        outgoing[from] += route.mass;
        incoming[to] += route.mass;
        downstream[from].push(to);
    }
    while let Some(from) = pending.pop() {
        poll(gate)?;
        for &to in &downstream[from] {
            if !reachable[to] {
                reachable[to] = true;
                pending.push(to);
            }
        }
    }
    if reachable.iter().any(|yes| !yes) {
        return Err(Error::Invalid(
            "mass-route component has no initial enthalpy storage",
        ));
    }
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        let supplied = old[i] + incoming[i];
        let removed = new[i] + outgoing[i];
        if !supplied.is_finite()
            || !removed.is_finite()
            || removed <= 0.0
            || (supplied - removed).abs() > 2048.0 * f64::EPSILON * supplied.max(removed)
            || !b[i].is_finite()
        {
            return Err(Error::Invalid(
                "native cell mass routes do not close to roundoff",
            ));
        }
        coo.push(i, i, removed);
    }
    for route in &trace.routes {
        coo.push(index[&route.to], index[&route.from], -route.mass);
    }
    let a = crate::conjugate::scale_rows(&coo, &mut b);
    if b.iter().any(|x| !x.is_finite()) {
        return Err(Error::Invalid("scaled enthalpy transport overflows"));
    }
    let outcome = bicgstab_ilu0(
        "free-surface enthalpy",
        &a,
        &b,
        &mut h,
        controls.relative_tolerance,
        controls.max_iterations,
        gate,
    )?;
    if h.iter().any(|value| !value.is_finite()) || !outcome.relative_residual.is_finite() {
        return Err(Error::Invalid("nonfinite accepted enthalpy transport"));
    }
    let result = index.into_iter().map(|(key, i)| (key, h[i])).collect();
    Ok((result, outcome))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::freesurface3::transport::MassRoute3;

    #[test]
    fn g1_resolved_native_fills_advect_their_real_correction_routes_without_heat_creation() {
        let (mut lattice, mut trace, a, b) = FreeSurface3::fill_projection_fixture(1000);
        let before = trace.before.clone();
        let old_h: BTreeMap<_, _> = before
            .keys()
            .filter_map(|&location| {
                let p = FreeSurface3::transport_coordinates(location)?;
                Some((
                    location,
                    if p[0] >= 14 {
                        100.0
                    } else if location == b {
                        20.0
                    } else {
                        10.0
                    },
                ))
            })
            .collect();
        let energy_before = total(&before, &old_h, 1.0).unwrap();
        let gate = CancelGate::new_clock_free();
        lattice.resolve_interface_fills(&mut trace, &gate).unwrap();
        assert_eq!(trace.before, before);
        assert_eq!(trace.routes.len(), 2);
        let (h, outcome) = advect(
            &trace,
            &old_h,
            FreeSurfaceEnthalpyControls::default(),
            &gate,
        )
        .unwrap();
        assert!(outcome.relative_residual < 1e-11);
        assert!((h[&a] - 15.0).abs() < 1e-10);
        let returned = trace.routes[1].mass;
        let expected_b = (before[&b] * 20.0 + returned * 15.0) / trace.routes[0].mass;
        assert!((h[&b] - expected_b).abs() < 1e-10);
        assert_eq!(trace.after[&b], 0.0);
        for (&location, &value) in &h {
            if FreeSurface3::transport_coordinates(location).is_some_and(|p| p[0] >= 14) {
                assert!((value - 100.0).abs() < 1e-10);
                assert_eq!(
                    trace.after[&location].to_bits(),
                    before[&location].to_bits()
                );
            }
        }
        assert!((total(&trace.after, &h, 1.0).unwrap() - energy_before).abs() < 1e-9);
    }

    #[test]
    fn g1_local_conversion_routes_preserve_disconnected_drop_enthalpy() {
        let a = MassLocation3::Cell(0, 0);
        let b = MassLocation3::Cell(0, 1);
        let c = MassLocation3::Cell(64, 0);
        let d = MassLocation3::Cell(64, 1);
        let carry = MassLocation3::Carry;
        let gate = CancelGate::new_clock_free();
        let controls = FreeSurfaceEnthalpyControls::default();
        // Two drops have independent temperatures and simultaneous conversion
        // excess. A shared zero-storage pool would mix their outgoing enthalpy.
        // The unlocated old carry is real storage, but has no physical route.
        let mut trace = MassTransport3::new(
            BTreeMap::from([(a, 2.0), (b, 1.0), (c, 2.0), (d, 1.0), (carry, 0.5)]),
            8,
        );
        trace.after = BTreeMap::from([(a, 1.0), (b, 2.0), (c, 1.0), (d, 2.0), (carry, 0.5)]);
        trace.routes = vec![
            MassRoute3 {
                from: a,
                to: b,
                mass: 1.0,
            },
            MassRoute3 {
                from: c,
                to: d,
                mass: 1.0,
            },
        ];
        let old_h = BTreeMap::from([
            (a, 10.0),
            (b, 40.0),
            (c, 100.0),
            (d, 400.0),
            (carry, 1000.0),
        ]);
        let (h, _) = advect(&trace, &old_h, controls, &gate).unwrap();
        for (key, expected) in [
            (a, 10.0),
            (b, 25.0),
            (c, 100.0),
            (d, 250.0),
            (carry, 1000.0),
        ] {
            assert!((h[&key] - expected).abs() < 1e-10, "{key:?}: {}", h[&key]);
        }
        assert!((total(&trace.after, &h, 1.0).unwrap() - 1160.0).abs() < 1e-10);

        // A disappearing cell oversends mass and its adjacent recipient
        // supplies the conversion debt. Each local implicit cycle is anchored
        // by its own old storage, with no cross-drop equilibration.
        trace.before = BTreeMap::from([(a, 1.0), (b, 2.0), (c, 1.0), (d, 2.0), (carry, 0.5)]);
        trace.after = BTreeMap::from([(a, 0.0), (b, 3.0), (c, 0.0), (d, 3.0), (carry, 0.5)]);
        trace.routes = vec![
            MassRoute3 {
                from: a,
                to: b,
                mass: 2.0,
            },
            MassRoute3 {
                from: b,
                to: a,
                mass: 1.0,
            },
            MassRoute3 {
                from: c,
                to: d,
                mass: 2.0,
            },
            MassRoute3 {
                from: d,
                to: c,
                mass: 1.0,
            },
        ];
        let old_h = BTreeMap::from([
            (a, 40.0),
            (b, 10.0),
            (c, 400.0),
            (d, 100.0),
            (carry, 1000.0),
        ]);
        let (h, _) = advect(&trace, &old_h, controls, &gate).unwrap();
        for (key, expected) in [
            (a, 30.0),
            (b, 20.0),
            (c, 300.0),
            (d, 200.0),
            (carry, 1000.0),
        ] {
            assert!((h[&key] - expected).abs() < 1e-10, "{key:?}: {}", h[&key]);
        }
        assert!((total(&trace.after, &h, 1.0).unwrap() - 1160.0).abs() < 1e-10);
    }

    #[test]
    fn g1_exact_conversion_pool_and_reverse_debt_transport() {
        let a = MassLocation3::Cell(0, 0);
        let b = MassLocation3::Cell(0, 1);
        let pool = MassLocation3::Carry;
        let gate = CancelGate::new_clock_free();
        let controls = FreeSurfaceEnthalpyControls::default();
        // Positive conversion excess: one mass unit at h=10 passes through
        // the zero-storage pool into a recipient initially at h=40.
        let mut trace = MassTransport3::new(BTreeMap::from([(a, 2.0), (b, 1.0), (pool, 0.0)]), 3);
        trace.after = BTreeMap::from([(a, 1.0), (b, 2.0), (pool, 0.0)]);
        trace.routes = vec![
            MassRoute3 {
                from: a,
                to: pool,
                mass: 1.0,
            },
            MassRoute3 {
                from: pool,
                to: b,
                mass: 1.0,
            },
        ];
        let old_h = BTreeMap::from([(a, 10.0), (b, 40.0)]);
        let (h, _) = advect(&trace, &old_h, controls, &gate).unwrap();
        for (key, expected) in [(a, 10.0), (b, 25.0), (pool, 10.0)] {
            assert!((h[&key] - expected).abs() < 1e-10);
        }
        assert!((total(&trace.after, &h, 1.0).unwrap() - 60.0).abs() < 1e-10);

        // A disappearing interface sent two units but held only one. Native
        // negative conversion excess is pool->a, while negative redistribution
        // is b->pool. The implicit cycle has an exact donor state h_a=30 and
        // final recipient/pool h=20; the zero-mass donor is not retained heat.
        trace.before = BTreeMap::from([(a, 1.0), (b, 2.0), (pool, 0.0)]);
        trace.after = BTreeMap::from([(a, 0.0), (b, 3.0), (pool, 0.0)]);
        trace.routes = vec![
            MassRoute3 {
                from: a,
                to: b,
                mass: 2.0,
            },
            MassRoute3 {
                from: pool,
                to: a,
                mass: 1.0,
            },
            MassRoute3 {
                from: b,
                to: pool,
                mass: 1.0,
            },
        ];
        let old_h = BTreeMap::from([(a, 40.0), (b, 10.0)]);
        let (h, _) = advect(&trace, &old_h, controls, &gate).unwrap();
        for (key, expected) in [(a, 30.0), (b, 20.0), (pool, 20.0)] {
            assert!((h[&key] - expected).abs() < 1e-10);
        }
        assert!((total(&trace.after, &h, 1.0).unwrap() - 60.0).abs() < 1e-10);
    }
}
