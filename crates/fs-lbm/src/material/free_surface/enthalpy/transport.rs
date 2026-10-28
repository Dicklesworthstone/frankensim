//! Conservative implicit donor transport over native transferred masses.
//!
//! Row i is `(m_new_i + sum(outgoing dm))*h_i - sum(incoming dm*h_donor)
//! = m_old_i*h_old_i`. The original mass balance makes this an anchored
//! M-matrix when each participating component receives old positive storage.
//! Conversion carry is an ordinary node with its original old/new mass. No
//! artificial mass, heat, concentration floor or diagonal regularizer is used.
//!
//! Solve for h minus one old-storage reference per connected route component.
//! Otherwise an arbitrary enthalpy datum enters the Krylov relative tolerance
//! and can suppress actual temperature mixing. The native mass-closure defect
//! remains explicitly in the shifted right-hand side; centering never repairs
//! a route or silently discards its associated reference-energy contribution.

use super::*;
use crate::conjugate::krylov::{bicgstab_ilu0, residual_norms};
use crate::d3q19::freesurface3::transport::MassTransport3;
use fs_sparse::Coo;

#[derive(Clone, Copy, Debug)]
pub(super) struct TransportOutcome {
    pub iterations: usize,
    pub relative_residual: f64,
    pub absolute_residual_j_kg: f64,
}

fn component_root(parents: &mut [usize], mut i: usize) -> usize {
    while parents[i] != i {
        parents[i] = parents[parents[i]];
        i = parents[i];
    }
    i
}

pub(super) fn advect(
    trace: &MassTransport3,
    old_h: &BTreeMap<MassLocation3, f64>,
    controls: FreeSurfaceEnthalpyControls,
    gate: &CancelGate,
) -> Result<(BTreeMap<MassLocation3, f64>, TransportOutcome), Error> {
    if !controls.absolute_transport_tolerance_j_kg.is_finite()
        || controls.absolute_transport_tolerance_j_kg < 0.0
    {
        return Err(Error::Invalid(
            "invalid absolute enthalpy transport tolerance",
        ));
    }
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
    let mut components: Vec<_> = (0..n).collect();
    let mut component_sizes = vec![1usize; n];
    for (&key, &i) in &index {
        poll(gate)?;
        old[i] = trace.before.get(&key).copied().unwrap_or(0.0);
        new[i] = trace.after.get(&key).copied().unwrap_or(0.0);
        if old[i] > 0.0 {
            h[i] = *old_h.get(&key).ok_or(Error::Invalid(
                "positive source mass needs retained enthalpy",
            ))?;
            if !h[i].is_finite() {
                return Err(Error::Invalid("nonfinite source enthalpy"));
            }
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
        let mut a = component_root(&mut components, from);
        let mut b = component_root(&mut components, to);
        if a != b {
            if component_sizes[a] < component_sizes[b] {
                core::mem::swap(&mut a, &mut b);
            }
            components[b] = a;
            component_sizes[a] += component_sizes[b];
        }
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
    let mut references = vec![None; n];
    for i in 0..n {
        poll(gate)?;
        components[i] = component_root(&mut components, i);
        if old[i] > 0.0 {
            // Index order is the deterministic native-location order. A born
            // cell or empty carry must inherit the component's real storage,
            // not introduce a zero-enthalpy reservoir into the linear solve.
            references[components[i]].get_or_insert(h[i]);
        }
    }
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        let supplied = old[i] + incoming[i];
        let removed = new[i] + outgoing[i];
        if !supplied.is_finite()
            || !removed.is_finite()
            || removed <= 0.0
            || (supplied - removed).abs() > 2048.0 * f64::EPSILON * supplied.max(removed)
        {
            return Err(Error::Invalid(
                "native cell mass routes do not close to roundoff",
            ));
        }
        let reference = references[components[i]].expect("reachable old storage");
        h[i] = if old[i] > 0.0 { h[i] - reference } else { 0.0 };
        // A(h_ref + delta_h) = m_old*h_old. In exact continuity the
        // second term vanishes. Retaining it preserves the original energy
        // equation even for the native, admitted roundoff mass discrepancy.
        b[i] = (supplied - removed).mul_add(reference, old[i] * h[i]);
        if !h[i].is_finite() || !b[i].is_finite() {
            return Err(Error::Invalid("relative enthalpy transport overflows"));
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
    for i in 0..n {
        h[i] += references[components[i]].expect("reachable old storage");
    }
    if h.iter().any(|value| !value.is_finite()) || !outcome.relative_residual.is_finite() {
        return Err(Error::Invalid("nonfinite accepted enthalpy transport"));
    }
    // Restoring a large datum can erase a correctly solved physical change.
    // Audit the enthalpies that will actually be stored, not the unrounded
    // centered iterate. The absolute physical allowance is explicit; an
    // arbitrary datum never supplies an implicit error scale.
    let represented: Vec<_> = h
        .iter()
        .enumerate()
        .map(|(i, &value)| value - references[components[i]].expect("reachable old storage"))
        .collect();
    poll(gate)?;
    let (absolute_residual_j_kg, rhs_norm) = residual_norms(&a, &b, &represented);
    poll(gate)?;
    let relative_residual = if rhs_norm == 0.0 && absolute_residual_j_kg == 0.0 {
        0.0
    } else {
        absolute_residual_j_kg / rhs_norm
    };
    let allowance =
        controls.absolute_transport_tolerance_j_kg + controls.relative_tolerance * rhs_norm;
    if !absolute_residual_j_kg.is_finite()
        || !rhs_norm.is_finite()
        || !relative_residual.is_finite()
        || !allowance.is_finite()
        || absolute_residual_j_kg > allowance
    {
        return Err(crate::conjugate::ChtError::SolverNotConverged {
            system: "free-surface represented enthalpy",
            iterations: outcome.iterations,
            relative_residual,
            tolerance: if rhs_norm > 0.0 {
                allowance / rhs_norm
            } else {
                controls.relative_tolerance
            },
        }
        .into());
    }
    let result = index.into_iter().map(|(key, i)| (key, h[i])).collect();
    Ok((
        result,
        TransportOutcome {
            iterations: outcome.iterations,
            relative_residual,
            absolute_residual_j_kg,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::d3q19::freesurface3::transport::MassRoute3;

    fn donor_pair() -> (MassTransport3, MassLocation3, MassLocation3) {
        let a = MassLocation3::Cell(0, 0);
        let b = MassLocation3::Cell(0, 1);
        let mut trace = MassTransport3::new(BTreeMap::from([(a, 2.0), (b, 1.0)]), 1);
        trace.after = BTreeMap::from([(a, 1.0), (b, 2.0)]);
        trace.routes.push(MassRoute3 {
            from: a,
            to: b,
            mass: 1.0,
        });
        (trace, a, b)
    }

    #[test]
    fn g1_donor_mixing_transports_latent_reference_once_at_every_energy_datum() {
        let (trace, a, b) = donor_pair();
        let gate = CancelGate::new_clock_free();
        // A synthetic liquid with Cp=4096 J/(kg K), carrying the already
        // absorbed fusion enthalpy. No new latent source acts during mixing.
        let fusion = 334_000.0;
        for datum in [0.0, 2.0f64.powi(50), -2.0f64.powi(50)] {
            let reference = datum + fusion;
            let old = BTreeMap::from([(a, reference + 4096.0), (b, reference + 12_288.0)]);
            let (h, outcome) =
                advect(&trace, &old, FreeSurfaceEnthalpyControls::default(), &gate).unwrap();
            // One unit at 274.15 K joins one unit at 276.15 K. The remaining
            // donor stays at 274.15 K; the recipient must reach 275.15 K.
            assert_eq!(h[&a] - reference, 4096.0);
            assert_eq!(h[&b] - reference, 8192.0);
            let recipient_temperature = 273.15 + (h[&b] - reference) / 4096.0;
            assert_eq!(recipient_temperature, 275.15);
            assert!(
                outcome.iterations > 0,
                "a large reference must not skip mixing"
            );
            assert!(
                outcome.relative_residual
                    <= FreeSurfaceEnthalpyControls::default().relative_tolerance
            );
            assert_eq!(
                trace.after[&a] * (h[&a] - reference) + trace.after[&b] * (h[&b] - reference),
                trace.before[&a] * (old[&a] - reference) + trace.before[&b] * (old[&b] - reference),
            );
        }
    }

    #[test]
    fn g3_disconnected_storage_cannot_hide_another_drops_mixing_residual() {
        let (mut trace, a, b) = donor_pair();
        let c = MassLocation3::Cell(64, 0);
        let d = MassLocation3::Cell(64, 1);
        let carry = MassLocation3::Carry;
        trace.before.extend([(c, 1.0), (d, 3.0), (carry, 2.0)]);
        trace.after.extend([(c, 0.0), (d, 4.0), (carry, 2.0)]);
        trace.routes.push(MassRoute3 {
            from: c,
            to: d,
            mass: 1.0,
        });
        let reference = 2.0f64.powi(50);
        let old = BTreeMap::from([
            (a, 10.0),
            (b, 40.0),
            (c, reference + 100.0),
            (d, reference + 400.0),
            (carry, -reference),
        ]);
        let (h, _) = advect(
            &trace,
            &old,
            FreeSurfaceEnthalpyControls::default(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
        assert!((h[&a] - 10.0).abs() < 1e-10);
        assert!((h[&b] - 25.0).abs() < 1e-10);
        assert_eq!(h[&c] - reference, 100.0);
        assert_eq!(h[&d] - reference, 325.0);
        assert_eq!(h[&carry], -reference);
    }

    #[test]
    fn g1_newborn_cell_and_empty_conversion_carry_inherit_the_donor_reference() {
        let a = MassLocation3::Cell(0, 0);
        let b = MassLocation3::Cell(0, 1);
        let carry = MassLocation3::Carry;
        let mut trace = MassTransport3::new(BTreeMap::from([(a, 2.0)]), 2);
        trace.after = BTreeMap::from([(a, 0.0), (b, 2.0), (carry, 0.0)]);
        trace.routes = vec![
            MassRoute3 {
                from: a,
                to: carry,
                mass: 2.0,
            },
            MassRoute3 {
                from: carry,
                to: b,
                mass: 2.0,
            },
        ];
        for reference in [0.0, 2.0f64.powi(50), -2.0f64.powi(50)] {
            let old = BTreeMap::from([(a, reference + 4096.0)]);
            let (h, _) = advect(
                &trace,
                &old,
                FreeSurfaceEnthalpyControls::default(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
            for location in [a, b, carry] {
                assert_eq!(h[&location], old[&a]);
            }
        }
    }

    #[test]
    fn g0_unrepresentable_fractional_mixing_refuses_the_published_enthalpy() {
        let (mut trace, a, b) = donor_pair();
        trace.before.insert(b, 2.0);
        trace.after.insert(b, 3.0);
        let controls = FreeSurfaceEnthalpyControls::default();
        let gate = CancelGate::new_clock_free();
        let (h, outcome) = advect(
            &trace,
            &BTreeMap::from([(a, 1.0), (b, 0.0)]),
            controls,
            &gate,
        )
        .unwrap();
        assert!((h[&b] - 1.0 / 3.0).abs() < 1e-15);
        assert!(outcome.relative_residual <= controls.relative_tolerance);
        for reference in [2.0f64.powi(52), -2.0f64.powi(52)] {
            let old = BTreeMap::from([(a, reference + 1.0), (b, reference)]);
            let refusal = advect(&trace, &old, controls, &gate).unwrap_err();
            assert!(
                matches!(refusal, Error::Solve(crate::conjugate::ChtError::SolverNotConverged {
                system: "free-surface represented enthalpy", relative_residual, ..
            }) if relative_residual > controls.relative_tolerance)
            );
        }
    }

    #[test]
    fn g0_explicit_absolute_allowance_handles_resolved_native_roundoff() {
        let a = MassLocation3::Cell(0, 0);
        let mut trace = MassTransport3::new(BTreeMap::from([(a, 1.0)]), 0);
        trace.after = BTreeMap::from([(a, 1.0 + f64::EPSILON)]);
        let old = BTreeMap::from([(a, 1e6)]);
        let controls = FreeSurfaceEnthalpyControls::default();
        let gate = CancelGate::new_clock_free();
        let (h, outcome) = advect(&trace, &old, controls, &gate).unwrap();
        let exact_increment = -f64::EPSILON * 1e6 / (1.0 + f64::EPSILON);
        let expected_residual = ((h[&a] - old[&a]) - exact_increment).abs();
        assert!((outcome.absolute_residual_j_kg - expected_residual).abs() < 1e-24);
        assert!(outcome.absolute_residual_j_kg > 0.0);
        assert!(outcome.absolute_residual_j_kg <= controls.absolute_transport_tolerance_j_kg);
        assert!(outcome.relative_residual > controls.relative_tolerance);
        // The strict relative-only choice is still available explicitly.
        let refusal = advect(
            &trace,
            &old,
            FreeSurfaceEnthalpyControls {
                absolute_transport_tolerance_j_kg: 0.0,
                ..controls
            },
            &gate,
        )
        .unwrap_err();
        assert!(matches!(
            refusal,
            Error::Solve(crate::conjugate::ChtError::SolverNotConverged {
                system: "free-surface represented enthalpy",
                ..
            })
        ));
        for invalid in [-1.0, f64::INFINITY, f64::NAN] {
            assert!(matches!(
                advect(
                    &trace,
                    &old,
                    FreeSurfaceEnthalpyControls {
                        absolute_transport_tolerance_j_kg: invalid,
                        ..controls
                    },
                    &gate
                ),
                Err(Error::Invalid(_))
            ));
        }
    }

    #[test]
    fn g0_centering_preserves_native_mass_defects_instead_of_repairing_them() {
        let (mut trace, a, b) = donor_pair();
        let reference = 2.0f64.powi(40);
        let old = BTreeMap::from([(a, reference), (b, reference + 30.0)]);
        let gate = CancelGate::new_clock_free();
        let controls = FreeSurfaceEnthalpyControls::default();
        let mismatch = 2.0f64.powi(-42);
        *trace.after.get_mut(&b).unwrap() += mismatch;
        let (h, _) = advect(&trace, &old, controls, &gate).unwrap();
        // The explicitly retained mass defect contributes -0.25 energy units
        // in the centered donor equation; dropping it would incorrectly give 15.
        let expected_increment = (30.0 - mismatch * reference) / (2.0 + mismatch);
        assert!((h[&b] - reference - expected_increment).abs() < 0.00025);
        assert!((h[&b] - reference - 15.0).abs() > 0.1);
        *trace.after.get_mut(&b).unwrap() += 1e-8;
        assert!(advect(&trace, &old, controls, &gate).is_err());
    }

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
