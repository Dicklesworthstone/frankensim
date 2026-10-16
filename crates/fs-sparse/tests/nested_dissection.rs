//! IRP .67 G0/G1: actual LDLT owner, analytic spectra and measured fill.
use fs_sparse::{
    Coo, Csr, DirectOrdering, LdltError, LdltOptions, NdOptions, SeparatorTree, SymbolicLdlt,
    nd_order, nd_order_with_coordinates,
};
use std::collections::BTreeSet;
use std::time::Instant;

fn policy() -> NdOptions {
    NdOptions {
        leaf_size: 64,
        refinement_passes: 4,
    }
}
fn grid(s: usize, shift: f64) -> Csr {
    let n = s * s * s;
    let mut rp = Vec::with_capacity(n + 1);
    let mut ci = Vec::new();
    let mut va = Vec::new();
    rp.push(0);
    for z in 0..s {
        for y in 0..s {
            for x in 0..s {
                let v = (z * s + y) * s + x;
                let mut row = vec![(v, 6.0 - shift)];
                for (b, limit, stride) in [(x, s, 1), (y, s, s), (z, s, s * s)] {
                    if b > 0 {
                        row.push((v - stride, -1.0));
                    }
                    if b + 1 < limit {
                        row.push((v + stride, -1.0));
                    }
                }
                row.sort_by_key(|p| p.0);
                for (c, a) in row {
                    ci.push(c);
                    va.push(a);
                }
                rp.push(ci.len());
            }
        }
    }
    Csr::from_parts(n, n, rp, ci, va)
}
fn xyz(s: usize) -> Vec<[f64; 3]> {
    (0..s * s * s)
        .map(|v| [(v % s) as f64, ((v / s) % s) as f64, (v / (s * s)) as f64])
        .collect()
}
fn validate(a: &Csr, t: &SeparatorTree) -> BTreeSet<usize> {
    match t {
        SeparatorTree::Leaf(v) => v.iter().copied().collect(),
        SeparatorTree::Branch {
            left,
            right,
            separator,
        } => {
            let mut l = validate(a, left);
            let r = validate(a, right);
            let s: BTreeSet<_> = separator.iter().copied().collect();
            assert!(!l.is_empty() && !r.is_empty());
            assert!(l.is_disjoint(&r) && l.is_disjoint(&s) && r.is_disjoint(&s));
            for &v in &l {
                assert!(
                    a.row(v).0.iter().all(|u| !r.contains(u)),
                    "crossing edge survived separator"
                );
            }
            l.extend(r);
            l.extend(s);
            l
        }
    }
}
fn check(a: &Csr, coords: Option<&[[f64; 3]]>) {
    let ordering = nd_order_with_coordinates(a, &policy(), coords).unwrap();
    let mut sorted = ordering.permutation.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, (0..a.nrows()).collect::<Vec<_>>());
    assert_eq!(validate(a, &ordering.tree), (0..a.nrows()).collect());
    fn positions(t: &SeparatorTree, p: &[usize]) -> Vec<usize> {
        match t {
            SeparatorTree::Leaf(v) => v.iter().map(|&i| p[i]).collect(),
            SeparatorTree::Branch {
                left,
                right,
                separator,
            } => {
                let mut v = positions(left, p);
                v.extend(positions(right, p));
                for &s in separator {
                    assert!(
                        v.iter().all(|&i| i < p[s]),
                        "separator must be eliminated last"
                    );
                }
                v.extend(separator.iter().map(|&i| p[i]));
                v
            }
        }
    }
    let mut inverse = vec![0; a.nrows()];
    for (k, &v) in ordering.permutation.iter().enumerate() {
        inverse[v] = k;
    }
    positions(&ordering.tree, &inverse);
}
#[test]
fn g0_nd_is_a_permutation_and_separators_separate() {
    check(&grid(7, 0.0), None);
    check(&grid(7, 0.0), Some(&xyz(7)));
    let mut c = Coo::new(160, 160);
    for i in 0..160 {
        c.push(i, i, 4.0);
        for j in i + 1..160 {
            if (i * 79 + j * 31) % 137 < 3 {
                c.push(i, j, -1.0);
                c.push(j, i, -1.0);
            }
        }
    }
    check(&c.assemble(), None);
    check(&Csr::identity(0), None);
    check(&Csr::identity(129), None);
}
#[test]
fn g0_nd_deterministic_repeat_bitwise() {
    let a = grid(9, 0.0);
    for coords in [None, Some(xyz(9))] {
        let x = nd_order_with_coordinates(&a, &policy(), coords.as_deref()).unwrap();
        for _ in 0..3 {
            assert_eq!(
                x,
                nd_order_with_coordinates(&a, &policy(), coords.as_deref()).unwrap()
            );
        }
    }
}
#[test]
fn g0_nd_leaf_size_bounds() {
    let a = grid(3, 0.0);
    assert!(
        nd_order(
            &a,
            &NdOptions {
                leaf_size: 0,
                refinement_passes: 4
            }
        )
        .is_err()
    );
    for leaf_size in [1, a.nrows()] {
        let o = nd_order(
            &a,
            &NdOptions {
                leaf_size,
                refinement_passes: 4,
            },
        )
        .unwrap();
        assert_eq!(validate(&a, &o.tree).len(), a.nrows());
    }
    assert!(nd_order_with_coordinates(&a, &policy(), Some(&[[f64::NAN; 3]; 27])).is_err());
    assert!(nd_order_with_coordinates(&a, &policy(), Some(&[[0.0; 3]; 26])).is_err());
    for refinement_passes in [0, 64] {
        nd_order(
            &a,
            &NdOptions {
                leaf_size: 1,
                refinement_passes,
            },
        )
        .unwrap();
    }
    assert!(
        nd_order(
            &a,
            &NdOptions {
                leaf_size: 1,
                refinement_passes: 65
            }
        )
        .is_err()
    );
    let rectangular = Coo::new(2, 3).assemble();
    assert!(matches!(
        nd_order(&rectangular, &policy()),
        Err(LdltError::NotSquare { .. })
    ));
    let mut asymmetric = Coo::new(2, 2);
    asymmetric.push(0, 1, 1.0);
    assert!(matches!(
        nd_order(&asymmetric.assemble(), &policy()),
        Err(LdltError::StructurallyAsymmetric { .. })
    ));
}
#[test]
fn g0_factor_budget_cap_and_cap_plus_one() {
    let a = grid(4, 0.0);
    let sym = SymbolicLdlt::analyze(&a, DirectOrdering::NestedDissection(policy())).unwrap();
    let bytes = sym.estimated_factor_bytes();
    let options = LdltOptions {
        max_factor_bytes: Some(bytes),
        ..LdltOptions::default()
    };
    let f = sym.factor(&a, &options).unwrap();
    assert!(f.stats().peak_factor_bytes <= bytes);
    let b = LdltOptions {
        max_factor_bytes: Some(bytes - 1),
        ..options
    };
    assert!(
        matches!(sym.factor(&a,&b),Err(LdltError::FactorBudget {estimated_bytes,cap}) if estimated_bytes==bytes && cap+1==bytes)
    );
    // Budget refusal precedes numeric-value inspection, proving that the
    // cap+1 path did not merely factor and inspect its eventual statistics.
    let mut bad = Coo::new(1, 1);
    bad.push(0, 0, f64::NAN);
    let bad = bad.assemble();
    let s = SymbolicLdlt::analyze(&bad, DirectOrdering::Amd).unwrap();
    assert!(matches!(
        s.factor(
            &bad,
            &LdltOptions {
                max_factor_bytes: Some(0),
                ..options
            }
        ),
        Err(LdltError::FactorBudget { .. })
    ));
}
#[test]
fn g1_nd_inertia_matches_analytic_3d_grid() {
    let s = 7;
    let a = grid(s, 0.0);
    let coordinates = xyz(s);
    let sym = SymbolicLdlt::analyze_with_coordinates(
        &a,
        DirectOrdering::NestedDissection(policy()),
        Some(&coordinates),
    )
    .unwrap();
    for shift in [-0.1, 0.731, 2.371, 5.123, 12.1] {
        let negative = (1..=s)
            .flat_map(|i| (1..=s).flat_map(move |j| (1..=s).map(move |k| (i, j, k))))
            .filter(|&(i, j, k)| {
                [i, j, k]
                    .iter()
                    .map(|&v| 2.0 - 2.0 * (std::f64::consts::PI * v as f64 / (s + 1) as f64).cos())
                    .sum::<f64>()
                    < shift
            })
            .count();
        let f = sym
            .factor(
                &grid(s, shift),
                &LdltOptions {
                    max_factor_bytes: Some(512 * 1024 * 1024),
                    ..LdltOptions::default()
                },
            )
            .unwrap();
        assert_eq!(f.inertia().negative, negative);
        assert_eq!(f.inertia().positive, s * s * s - negative);
    }
}
/// Release-only scale lane, never an ignored or reduced substitute for the
/// specified 30^3/45^3/60^3 fill comparison. No numeric factors are needed to
/// count symbolic fill; their independently measured flops are NA here.
#[test]
#[cfg(not(debug_assertions))]
fn g1_nd_fill_beats_amd_on_3d_grids() {
    let mut fill = Vec::new();
    let start = Instant::now();
    for s in [30, 45, 60] {
        let a = grid(s, 0.0);
        let coordinates = xyz(s);
        let timer = Instant::now();
        let amd = SymbolicLdlt::analyze(&a, DirectOrdering::Amd).unwrap();
        let amd_s = timer.elapsed().as_secs_f64();
        let timer = Instant::now();
        let nd = SymbolicLdlt::analyze_with_coordinates(
            &a,
            DirectOrdering::NestedDissection(policy()),
            Some(&coordinates),
        )
        .unwrap();
        let nd_s = timer.elapsed().as_secs_f64();
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.5-W1\",\"test\":\"g1_nd_fill_beats_amd_on_3d_grids\",\"inputs\":{{\"fixture\":\"grid3d-{s}\",\"n\":{},\"ordering\":\"geometric-nd-vs-amd\",\"leaf_size\":64}},\"observed\":{{\"amd_nnz_l\":{},\"nd_nnz_l\":{},\"estimated_flops_amd\":{},\"estimated_flops_nd\":{},\"numeric_peak_bytes\":null,\"symbolic_factor_payload_estimate\":{}}},\"receipt\":{{\"host_fingerprint\":\"{}-{}\",\"amd_wall_s\":{amd_s},\"nd_wall_s\":{nd_s}}},\"tolerance\":\"nnz(ND)<=nnz(AMD)\",\"verdict\":\"{}\",\"elapsed_s\":{}}}",
            a.nrows(),
            amd.nnz_l(),
            nd.nnz_l(),
            amd.estimated_flops(),
            nd.estimated_flops(),
            nd.estimated_factor_bytes(),
            std::env::consts::ARCH,
            std::env::consts::OS,
            if nd.nnz_l() <= amd.nnz_l() {
                "pass"
            } else {
                "fail"
            },
            start.elapsed().as_secs_f64()
        );
        assert!(
            nd.nnz_l() <= amd.nnz_l(),
            "{s}^3: ND{} > AMD{}",
            nd.nnz_l(),
            amd.nnz_l()
        );
        fill.push((a.nrows() as f64, nd.nnz_l() as f64));
    }
    let exponent = (fill[2].1 / fill[0].1).ln() / (fill[2].0 / fill[0].0).ln();
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.5-W1\",\"test\":\"nd-fill-growth\",\"observed\":{{\"exponent_vs_dof\":{exponent},\"reference_3d_exponent\":1.3333333333333333}},\"verdict\":\"reported\"}}"
    );
    assert!(exponent.is_finite());
}
