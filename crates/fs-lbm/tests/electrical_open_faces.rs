//! G0/G1/G3/G4 electrical topology through the original current/Joule solver.
//! Open faces are not high resistances and do not change bulk material identity.
use std::cell::Cell;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricalContact, ElectricSetup, Electrode,
    ElectrodeDrive, VoxelDomain};

fn bar(n: usize, drive: ElectrodeDrive) -> (VoxelDomain, ElectricSetup) {
    (VoxelDomain::new(n, 1, 1, 1.0).unwrap(), ElectricSetup {
        conductors: vec![], electrodes: vec![
            Electrode { lo: [0,0,0], hi: [1,1,1], drive },
            Electrode { lo: [n-1,0,0], hi: [n,1,1], drive: ElectrodeDrive::Voltage(0.0) },
        ],
    })
}
fn contact(cells: [usize; 2]) -> ElectricalContact {
    ElectricalContact { cells, resistance_ohm_m2: 3.0, heat_fraction_to_first: 0.25 }
}
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-9, "{a} != {b}"); }

#[test]
fn empty_openings_preserve_the_resistive_solution_exactly() {
    let (d, e) = bar(4, ElectrodeDrive::Voltage(9.0));
    let gate = CancelGate::new(); let faces = [contact([1,2])];
    let old = e.solve_with_cell_contacts(&d, &[300.0;4], &faces,
        |c,_,_| Ok(Some(if c==1 {2.0} else {4.0})), &gate).unwrap();
    let new = e.solve_with_cell_contact_topology(&d, &[300.0;4], &faces, &[],
        |c,_,_| Ok(Some(if c==1 {2.0} else {4.0})), &gate).unwrap();
    assert_eq!(old, new);
    near(new.electrodes[0].current_a, 1.0);
    near(new.total_joule_w, 9.0);
}

#[test]
fn opening_and_reclosing_only_changes_the_selected_electrical_face() {
    let (d, e) = bar(5, ElectrodeDrive::Voltage(10.0));
    let gate = CancelGate::new(); let faces = [contact([1,2])];
    let closed = e.solve_with_cell_contact_topology(&d, &[300.0;5], &faces, &[],
        |_,_,_| Ok(Some(2.0)), &gate).unwrap();
    near(closed.electrodes[0].current_a, 10.0/9.0);
    for opening in [[1,2], [2,1]] {
        let open = e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &[opening],
            |_,_,_| Ok(Some(2.0)), &gate).unwrap();
        for e in &open.electrodes { near(e.current_a, 0.0); }
        assert!(open.joule_w.iter().all(|q| q.abs() < 1e-25));
        near(open.total_joule_w, 0.0); near(open.delivered_w, 0.0);
        near(open.potential_v[1], 10.0); near(open.potential_v[2], 0.0);
    }
    let again = e.solve_with_cell_contact_topology(&d, &[300.0;5], &faces, &[],
        |_,_,_| Ok(Some(2.0)), &gate).unwrap();
    assert_eq!(closed, again);
    assert_eq!(faces[0], contact([1,2]));
}

#[test]
fn floating_bulk_is_still_evaluated_and_no_leakage_or_joule_heat_is_fabricated() {
    let (d, e) = bar(5, ElectrodeDrive::Voltage(10.0)); let calls = Cell::new(0);
    let gate = CancelGate::new(); let openings = [[0,1], [3,4]];
    let result = e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &openings,
        |_,_,_| { calls.set(calls.get()+1); Ok(Some(2.0)) }, &gate).unwrap();
    assert_eq!(calls.get(), 3);
    assert!(result.potential_v[1..4].iter().all(|v| v.is_nan()));
    assert_eq!(result.joule_w, vec![0.0;5]);
    assert_eq!(result.total_joule_w, 0.0); assert_eq!(result.delivered_w, 0.0);
    let failure = ChtError::InvalidInput { field: "test.material", reason: "unsupported floating material".into() };
    let attempt = e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &openings,
        |c,_,_| if c==2 { Err(failure.clone()) } else { Ok(Some(2.0)) }, &gate);
    assert_eq!(attempt.unwrap_err(), failure);
}

#[test]
fn a_disconnected_nonzero_current_source_refuses_but_zero_current_has_zero_work() {
    let (d, mut e) = bar(3, ElectrodeDrive::Current(1.0)); let gate = CancelGate::new();
    let result = e.solve_with_cell_contact_topology(&d, &[300.0;3], &[], &[[0,1]],
        |_,_,_| Ok(Some(2.0)), &gate);
    assert!(matches!(result, Err(ChtError::InvalidInput {field: "electric.electrodes", ..})));
    e.electrodes[0].drive = ElectrodeDrive::Current(0.0);
    let result = e.solve_with_cell_contact_topology(&d, &[300.0;3], &[], &[[0,1]],
        |_,_,_| Ok(Some(2.0)), &gate).unwrap();
    assert!(result.electrodes[0].voltage_v.is_nan());
    assert_eq!(result.electrodes[0].current_a, 0.0);
    assert_eq!(result.delivered_w, 0.0); assert_eq!(result.total_joule_w, 0.0);
}

#[test]
fn open_face_preserves_alternate_paths_and_matches_independent_kirchhoff_solution() {
    let d = VoxelDomain::new(3,2,1,1.0).unwrap(); let gate = CancelGate::new();
    for sign in [1.0, -1.0] {
        let e = ElectricSetup { conductors: vec![], electrodes: vec![
            Electrode {lo:[0,0,0],hi:[1,2,1],drive:ElectrodeDrive::Voltage(sign*8.0)},
            Electrode {lo:[2,0,0],hi:[3,2,1],drive:ElectrodeDrive::Voltage(0.0)},
        ]};
        let result = e.solve_with_cell_contact_topology(&d, &[300.0;6], &[], &[[0,1]],
            |_,_,_| Ok(Some(2.0)), &gate).unwrap();
        // Centre potentials obey 1.5 A - .5 B = 0; -.5 A + 2.5 B = 8.
        near(result.potential_v[1], sign*8.0/7.0);
        near(result.potential_v[4], sign*24.0/7.0);
        near(result.electrodes[0].current_a, sign*32.0/7.0);
        near(result.electrodes[1].current_a, -sign*32.0/7.0);
        near(result.total_joule_w, 256.0/7.0);
        near(result.delivered_w, result.total_joule_w);
        assert!(result.joule_w.iter().all(|q| *q >= 0.0));
    }
}

#[test]
fn conflicting_invalid_or_cancelled_topology_does_not_change_any_input() {
    let (d, e) = bar(5, ElectrodeDrive::Voltage(10.0)); let gate = CancelGate::new();
    for openings in [vec![[1,1]], vec![[0,2]], vec![[1,99]], vec![[1,2],[2,1]]] {
        assert!(e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &openings,
            |_,_,_| Ok(Some(2.0)), &gate).is_err());
    }
    assert!(e.solve_with_cell_contact_topology(&d, &[300.0;5], &[contact([1,2])], &[[2,1]],
        |_,_,_| Ok(Some(2.0)), &gate).is_err());
    let cancelled = CancelGate::new(); cancelled.request();
    assert!(matches!(e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &[[1,2]],
        |_,_,_| panic!("cancelled before material evaluation"), &cancelled), Err(ChtError::Cancelled)));
    let result = e.solve_with_cell_contact_topology(&d, &[300.0;5], &[], &[],
        |_,_,_| Ok(Some(2.0)), &gate).unwrap();
    near(result.electrodes[0].current_a, 10.0/6.0);
}
