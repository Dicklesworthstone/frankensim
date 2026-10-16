//! Battery for the fidelity-ladder registry (addendum Proposal 3). Covers
//! rung ordering + resolution, adjacency at the ends, the G0 transfer
//! approximation property (restrict∘prolongate = identity; prolongate∘restrict
//! idempotent), structured boundary errors (prolongate-at-top, restrict-at-
//! bottom, out-of-range rung, unknown kernel), the CHT ladder, and G5
//! determinism.

use fs_ladder::{AdjacentRungs, Ladder, LadderError, LadderRegistry, Refine1d, Transfer};

/// Bit-exact view of a float slice (avoids float-cmp lints and asserts
/// determinism to the bit).
fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

fn three_rung() -> Ladder {
    Ladder::new("demo", "coarse", 1.0, "bottom")
        .then(Box::new(Refine1d), "mid", 10.0, "middle")
        .then(Box::new(Refine1d), "fine", 100.0, "top")
}

#[test]
fn rungs_are_totally_ordered() {
    let l = three_rung();
    assert_eq!(l.len(), 3);
    assert_eq!(l.rung(0).unwrap().index, 0);
    assert_eq!(l.rung(1).unwrap().name, "mid");
    assert_eq!(l.rung(2).unwrap().index, 2);
    assert_eq!(l.bottom().name, "coarse");
    assert_eq!(l.top().name, "fine");
    // strictly increasing indices.
    for k in 0..l.len() {
        assert_eq!(l.rung(k).unwrap().index, k);
    }
}

#[test]
fn out_of_range_rung_is_a_structured_error() {
    let l = three_rung();
    match l.rung(9) {
        Err(LadderError::NoSuchRung { kernel, index, len }) => {
            assert_eq!(kernel, "demo");
            assert_eq!(index, 9);
            assert_eq!(len, 3);
        }
        other => panic!("expected NoSuchRung, got {other:?}"),
    }
    // the error teaches.
    assert!(l.rung(9).unwrap_err().to_string().contains("no rung 9"));
}

#[test]
fn adjacency_is_empty_in_the_right_direction_at_the_ends() {
    let l = three_rung();
    // bottom: no coarser, has finer.
    assert_eq!(
        l.adjacent_rungs(0).unwrap(),
        AdjacentRungs {
            coarser: None,
            finer: Some(1)
        }
    );
    // middle: both.
    assert_eq!(
        l.adjacent_rungs(1).unwrap(),
        AdjacentRungs {
            coarser: Some(0),
            finer: Some(2)
        }
    );
    // top: has coarser, no finer.
    assert_eq!(
        l.adjacent_rungs(2).unwrap(),
        AdjacentRungs {
            coarser: Some(1),
            finer: None
        }
    );
    // out of range is a structured error.
    assert!(matches!(
        l.adjacent_rungs(3),
        Err(LadderError::NoSuchRung { .. })
    ));
}

#[test]
fn transfer_g0_property_restrict_after_prolongate_is_identity() {
    let l = three_rung();
    let coarse = vec![0.0, 2.0, 4.0, 6.0];
    let fine = l.prolongate(0, &coarse).unwrap();
    // linear interpolation doubles-minus-one the sample count.
    assert_eq!(fine.len(), coarse.len() * 2 - 1);
    let back = l.restrict(1, &fine).unwrap();
    assert_eq!(
        bits(&back),
        bits(&coarse),
        "restrict∘prolongate must be identity"
    );
}

#[test]
fn transfer_prolongate_after_restrict_is_idempotent() {
    let t = Refine1d;
    let fine = vec![0.0, 5.0, 2.0, 9.0, 4.0];
    let pr = t.prolongate(&t.restrict(&fine));
    let pr2 = t.prolongate(&t.restrict(&pr));
    assert_eq!(
        bits(&pr2),
        bits(&pr),
        "prolongate∘restrict must be idempotent (a projection)"
    );
}

#[test]
fn prolongate_at_the_top_is_a_structured_error() {
    let l = three_rung();
    match l.prolongate(2, &[1.0, 2.0]) {
        Err(LadderError::AtTop { kernel, index }) => {
            assert_eq!(kernel, "demo");
            assert_eq!(index, 2);
        }
        other => panic!("expected AtTop, got {other:?}"),
    }
}

#[test]
fn restrict_at_the_bottom_is_a_structured_error() {
    let l = three_rung();
    match l.restrict(0, &[1.0, 2.0, 3.0]) {
        Err(LadderError::AtBottom { kernel, index }) => {
            assert_eq!(kernel, "demo");
            assert_eq!(index, 0);
        }
        other => panic!("expected AtBottom, got {other:?}"),
    }
}

#[test]
fn out_of_range_transfer_is_no_such_rung() {
    let l = three_rung();
    assert!(matches!(
        l.prolongate(5, &[1.0]),
        Err(LadderError::NoSuchRung { .. })
    ));
    assert!(matches!(
        l.restrict(5, &[1.0]),
        Err(LadderError::NoSuchRung { .. })
    ));
}

#[test]
fn registry_registers_and_resolves_kernels() {
    let mut r = LadderRegistry::new();
    r.register(three_rung());
    assert!(r.ladder("demo").is_ok());
    assert_eq!(r.kernels(), vec!["demo"]);
    // unknown kernel is a structured error.
    match r.ladder("ghost") {
        Err(LadderError::NoKernel { kernel }) => assert_eq!(kernel, "ghost"),
        other => panic!("expected NoKernel, got {other:?}"),
    }
}

#[test]
fn cht_ladder_makes_it_real() {
    // the correlation-based bottom rung from Proposal 7.
    let r = LadderRegistry::cht();
    let cht = r.ladder("cht").unwrap();
    assert_eq!(cht.len(), 3);
    assert_eq!(cht.bottom().name, "correlation-Nu");
    assert_eq!(cht.rung(1).unwrap().name, "RANS");
    assert_eq!(cht.top().name, "LES");
    // costs increase up the ladder (advisory).
    assert!(cht.bottom().relative_cost < cht.top().relative_cost);
    // the ladder actually transfers: prolongate then restrict recovers.
    let coarse = vec![1.0, 3.0, 5.0];
    let fine = cht.prolongate(0, &coarse).unwrap();
    assert_eq!(bits(&cht.restrict(1, &fine).unwrap()), bits(&coarse));
}

#[test]
fn refine1d_handles_degenerate_lengths() {
    let t = Refine1d;
    assert!(t.prolongate(&[]).is_empty());
    assert!(t.restrict(&[]).is_empty());
    // a single sample prolongates to itself and restricts to itself.
    assert_eq!(bits(&t.prolongate(&[7.0])), bits(&[7.0]));
    assert_eq!(bits(&t.restrict(&[7.0])), bits(&[7.0]));
}

#[test]
fn transfers_are_deterministic() {
    let l = three_rung();
    let coarse = vec![0.5, 1.5, 2.5];
    let a = l.prolongate(0, &coarse).unwrap();
    let b = l.prolongate(0, &coarse).unwrap();
    assert_eq!(
        bits(&a),
        bits(&b),
        "prolongate must be bit-identical on replay"
    );
    let ra = l.restrict(1, &a).unwrap();
    let rb = l.restrict(1, &b).unwrap();
    assert_eq!(bits(&ra), bits(&rb));
}
