use super::*;

fn shifted(x: [f64; 3], n: [f64; 3], h: f64) -> [f64; 3] {
    [x[0] + h * n[0], x[1] + h * n[1], x[2] + h * n[2]]
}

fn close(got: C64, want: C64, tol: f64) {
    assert!(
        (got - want).abs() <= tol * want.abs().max(1.0e-12),
        "{got:?} vs {want:?}"
    );
}

#[test]
fn helmholtz_point_kernels_g1_sign_gradient_and_mixed_derivative() {
    let x = [1.2, -0.4, 0.7];
    let y = [-0.2, 0.5, 0.1];
    let nx = [0.6, 0.0, 0.8];
    let ny = [0.0, 0.8, -0.6];
    for k in [0.0, 0.7, 4.0] {
        let single = HelmholtzSingle { k };
        let double = HelmholtzDouble { k };
        let adjoint = HelmholtzAdjointDouble { k };
        let mixed = HelmholtzHypersingular { k };
        let h = 1.0e-5;
        let fd_y =
            (single.eval(x, shifted(y, ny, h)) - single.eval(x, shifted(y, ny, -h))).scale(0.5 / h);
        let fd_x =
            (single.eval(shifted(x, nx, h), y) - single.eval(shifted(x, nx, -h), y)).scale(0.5 / h);
        close(double.eval_with_normal(x, y, ny), fd_y, 2.0e-8);
        close(adjoint.eval_with_normal(x, y, nx), fd_x, 2.0e-8);
        let fd_xy = (double.eval_with_normal(shifted(x, nx, h), y, ny)
            - double.eval_with_normal(shifted(x, nx, -h), y, ny))
        .scale(0.5 / h);
        close(mixed.eval_with_normals(x, y, nx, ny), fd_xy, 2.0e-8);
        close(
            double.eval_with_normal(x, y, ny),
            adjoint.eval_with_normal(y, x, ny),
            1.0e-14,
        );
        let mut normal_gradient = C64::ZERO;
        for a in 0..3 {
            let mut unit = [0.0; 3];
            unit[a] = 1.0;
            let component = HelmholtzGradient::new(k, a).unwrap().eval(x, y);
            let fd = (single.eval(shifted(x, unit, h), y) - single.eval(shifted(x, unit, -h), y))
                .scale(0.5 / h);
            close(component, fd, 2.0e-8);
            normal_gradient = normal_gradient + component.scale(nx[a]);
        }
        close(normal_gradient, adjoint.eval_with_normal(x, y, nx), 1.0e-14);
        assert_eq!(single.eval(x, x), C64::ZERO);
        assert_eq!(double.eval_with_normal(x, x, ny), C64::ZERO);
        assert_eq!(adjoint.eval_with_normal(x, x, nx), C64::ZERO);
        assert_eq!(mixed.eval_with_normals(x, x, nx, ny), C64::ZERO);
    }
    // Analytic outgoing phase at kr=pi/2, independently of the derivative code.
    let g = HelmholtzSingle {
        k: core::f64::consts::FRAC_PI_2,
    }
    .eval([1.0, 0.0, 0.0], [0.0; 3]);
    close(
        g,
        C64::new(0.0, 1.0 / (4.0 * core::f64::consts::PI)),
        1.0e-14,
    );
    assert!(HelmholtzGradient::new(1.0, 3).is_err());
    assert!(HelmholtzGradient::new(f64::NAN, 0).is_err());
    assert!(HelmholtzGradient::new(-1.0, 0).is_err());
}

#[test]
fn complex_tree_pair_partition_g0_no_missing_or_duplicate_interactions() {
    struct Constant;
    impl ComplexKernel for Constant {
        fn eval(&self, _: [f64; 3], _: [f64; 3]) -> C64 {
            C64::new(2.0, -3.0)
        }
    }
    // Three levels including parent M2L and separated leaves, so both M2M
    // and L2L influence the result (a depth-2 tree would not exercise this).
    let points = (0..72)
        .map(|i| {
            [
                (i % 9) as f64 / 8.0,
                ((i / 9) % 4) as f64 / 4.0,
                (i / 36) as f64 / 2.0,
            ]
        })
        .collect();
    let q: Vec<C64> = (0..72).map(|i| C64::new(i as f64 / 72.0, 0.25)).collect();
    let engine = FmmC64::new(&Constant, points, 3, 1).unwrap();
    assert!(!engine.interactions.is_empty());
    assert!(
        engine
            .interactions
            .iter()
            .any(|(a, _)| a.0 < engine.geometry.max_level)
    );
    assert!(
        engine
            .interactions
            .iter()
            .any(|(a, _)| a.0 == engine.geometry.max_level)
    );
    let got = engine.potentials(&q).unwrap();
    let direct = engine.direct(&q).unwrap();
    for (&g, &d) in got.iter().zip(&direct) {
        close(g, d, 2.0e-14);
    }
    let points = &engine.geometry.points;
    for t in 0..points.len() {
        for s in 0..points.len() {
            let leaf = |index| {
                engine
                    .geometry
                    .cells
                    .iter()
                    .find(|(_, c)| c.leaf && c.points.contains(&index))
                    .map(|(&k, _)| k)
                    .unwrap()
            };
            let (a, b) = (leaf(t), leaf(s));
            let mut count = usize::from(Fmm::adjacent(a, b));
            let ancestor =
                |k: CellKey, lv: u32| (lv, k.1 >> (k.0 - lv), k.2 >> (k.0 - lv), k.3 >> (k.0 - lv));
            for &(ma, mb) in &engine.interactions {
                if ma == ancestor(a, ma.0) && mb == ancestor(b, mb.0) {
                    count += 1;
                }
            }
            assert_eq!(count, 1, "pair ({t},{s}) has {count} owners");
        }
    }
}

#[test]
fn shared_m2l_blocks_match_per_pair_evaluation_g1() {
    // The same Helmholtz kernel without the invariance opt-in evaluates every
    // pair's block from its own grids; sharing may change only rounding.
    struct PerPair(HelmholtzSingle);
    impl ComplexKernel for PerPair {
        fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64 {
            self.0.eval(x, y)
        }
        fn wavenumber(&self) -> Option<f64> {
            self.0.wavenumber()
        }
    }
    let points: Vec<[f64; 3]> = (0..64)
        .map(|i| {
            let j = 0.1 * det::sin(i as f64);
            [
                (i % 4) as f64 + j,
                ((i / 4) % 4) as f64 - j,
                (i / 16) as f64 + 0.5 * j,
            ]
        })
        .collect();
    let q: Vec<C64> = (0..64)
        .map(|i| C64::new(det::cos(i as f64), 0.3 - i as f64 / 64.0))
        .collect();
    let offset = |(a, b): (CellKey, CellKey)| {
        let d = |i: u32, j: u32| i64::from(j) - i64::from(i);
        (a.0, d(a.1, b.1), d(a.2, b.2), d(a.3, b.3))
    };
    let single = HelmholtzSingle { k: 1.3 };
    let shared = FmmC64::new(&single, points.clone(), 3, 1).unwrap();
    let per_pair_kernel = PerPair(single);
    let opted_out = FmmC64::new(&per_pair_kernel, points.clone(), 3, 1).unwrap();
    assert_eq!(shared.interactions, opted_out.interactions);
    assert_eq!(opted_out.translations.len(), opted_out.interactions.len());
    // 64 occupied level-2 cells: 3096 pairs over at most 316 offsets.
    assert!(shared.translations.len() * 4 < shared.interactions.len());
    let mut owners = vec![0usize; shared.interactions.len()];
    for group in &shared.translations {
        for &index in group {
            assert_eq!(
                offset(shared.interactions[index]),
                offset(shared.interactions[group[0]])
            );
            owners[index] += 1;
        }
    }
    assert!(owners.iter().all(|&n| n == 1));
    // Rounding-level agreement, relative to the largest potential.
    let agree = |got: Vec<C64>, want: Vec<C64>| {
        let scale = want.iter().map(|w| w.abs()).fold(0.0, f64::max);
        for (g, w) in got.into_iter().zip(want) {
            assert!((g - w).abs() <= 1.0e-12 * scale, "{g:?} vs {w:?}");
        }
    };
    agree(
        shared.potentials(&q).unwrap(),
        opted_out.potentials(&q).unwrap(),
    );
    // Multi-channel operators index the shared row by channel; compare each
    // against the same engine forced back to one block per pair.
    let normals: Vec<[f64; 3]> = points
        .iter()
        .map(|p| {
            let v = [p[0] + 1.0, p[1] + 2.0, p[2] + 3.0];
            let n = det::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            [v[0] / n, v[1] / n, v[2] / n]
        })
        .collect();
    let build = |layer: usize| {
        let (p, n) = (points.clone(), normals.clone());
        match layer {
            0 => FmmC64::new_double(HelmholtzDouble { k: 1.3 }, p, n, 3, 1),
            1 => FmmC64::new_adjoint_double(HelmholtzAdjointDouble { k: 1.3 }, p, n, 3, 1),
            _ => FmmC64::new_hypersingular(HelmholtzHypersingular { k: 1.3 }, p, n, 3, 1),
        }
        .unwrap()
    };
    for layer in 0..3 {
        let shared = build(layer);
        let mut per_pair = build(layer);
        per_pair.translations = (0..per_pair.interactions.len()).map(|i| vec![i]).collect();
        assert!(shared.translations.len() * 4 < shared.interactions.len());
        agree(
            shared.potentials(&q).unwrap(),
            per_pair.potentials(&q).unwrap(),
        );
    }
}

#[test]
fn near_field_cache_is_bitwise_and_budgeted_g0() {
    let points: Vec<[f64; 3]> = (0..96)
        .map(|i| {
            let j = 0.07 * det::cos(i as f64);
            [
                (i % 6) as f64 + j,
                ((i / 6) % 4) as f64 - j,
                (i / 24) as f64,
            ]
        })
        .collect();
    let normals: Vec<[f64; 3]> = points
        .iter()
        .map(|p| {
            let v = [p[0] + 1.0, p[1] + 2.0, p[2] + 3.0];
            let n = det::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            [v[0] / n, v[1] / n, v[2] / n]
        })
        .collect();
    let q: Vec<C64> = (0..96)
        .map(|i| C64::new(det::sin(i as f64), 0.2 - i as f64 / 96.0))
        .collect();
    let build = || {
        FmmC64::new_hypersingular(
            HelmholtzHypersingular { k: 0.9 },
            points.clone(),
            normals.clone(),
            3,
            2,
        )
        .unwrap()
    };
    let plain = build();
    let mut cached = build();
    let bytes = cached.cache_near_field(usize::MAX).unwrap();
    assert!(bytes > 0 && bytes % core::mem::size_of::<C64>() == 0);
    let bits = |p: Vec<C64>| -> Vec<(u64, u64)> {
        p.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
    };
    assert_eq!(
        bits(plain.potentials(&q).unwrap()),
        bits(cached.potentials(&q).unwrap())
    );
    // The cap is inclusive; one byte less refuses and leaves the engine uncached.
    assert_eq!(build().cache_near_field(bytes).unwrap(), bytes);
    let mut short = build();
    assert!(matches!(
        short.cache_near_field(bytes - 1),
        Err(FmmError::WorkEnvelopeExceeded { requested, max, .. }) if requested == bytes && max == bytes - 1
    ));
    assert!(short.near.is_none());
    assert_eq!(
        bits(short.potentials(&q).unwrap()),
        bits(plain.potentials(&q).unwrap())
    );
}

#[test]
fn complex_invalid_requests_and_non_finite_kernels_g0() {
    let points = vec![[0.0; 3], [1.0, 0.0, 0.0]];
    let single = HelmholtzSingle { k: 1.0 };
    assert!(FmmC64::new(&single, Vec::new(), 4, 8).is_err());
    assert!(FmmC64::new(&single, vec![[f64::NAN; 3]], 4, 8).is_err());
    assert!(FmmC64::new(&single, points.clone(), 1, 8).is_err());
    assert!(FmmC64::new(&single, points.clone(), 13, 8).is_err());
    assert!(FmmC64::new(&single, points.clone(), 4, 0).is_err());
    for k in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(FmmC64::new(&HelmholtzSingle { k }, points.clone(), 4, 8).is_err());
    }
    for ns in [
        vec![],
        vec![[0.0; 3]; 2],
        vec![[f64::NAN; 3]; 2],
        vec![[2.0, 0.0, 0.0]; 2],
    ] {
        assert!(FmmC64::new_double(HelmholtzDouble { k: 1.0 }, points.clone(), ns, 4, 8).is_err());
    }
    let engine = FmmC64::new(&single, points.clone(), 4, 8).unwrap();
    assert!(engine.potentials(&[C64::ONE]).is_err());
    assert!(engine.direct(&[C64::ONE, C64::new(f64::NAN, 0.0)]).is_err());
    assert!(
        engine
            .potentials(&[C64::ONE, C64::new(0.0, f64::INFINITY)])
            .is_err()
    );
    struct Broken;
    impl ComplexKernel for Broken {
        fn eval(&self, _: [f64; 3], _: [f64; 3]) -> C64 {
            C64::new(f64::NAN, 0.0)
        }
    }
    let broken = FmmC64::new(&Broken, points, 4, 8).unwrap();
    assert!(matches!(
        broken.potentials(&[C64::ONE; 2]),
        Err(FmmError::NonFiniteOutput { .. })
    ));
    assert!(matches!(
        broken.direct(&[C64::ONE; 2]),
        Err(FmmError::NonFiniteOutput { .. })
    ));
    for (kd, tol) in [
        (f64::NAN, 1.0e-6),
        (-1.0, 1.0e-6),
        (1.0, 0.0),
        (1.0, f64::NAN),
        (1.0, 1.0),
    ] {
        assert!(helmholtz_order_for(kd, tol).is_err());
    }
}

#[test]
fn complex_normals_contract_exact_direct_kernel_g1() {
    let x = [0.3, 0.2, -0.1];
    let y = [-0.4, 0.7, 0.9];
    let nx = [0.6, 0.0, 0.8];
    let ny = [0.0, 0.8, -0.6];
    let points = vec![x, y];
    let normals = vec![nx, ny];
    let q = [C64::new(0.2, -0.7), C64::new(-0.4, 0.1)];
    let double = HelmholtzDouble { k: 2.0 };
    let adjoint = HelmholtzAdjointDouble { k: 2.0 };
    let mixed = HelmholtzHypersingular { k: 2.0 };
    let engines = [
        (
            FmmC64::new_double(double, points.clone(), normals.clone(), 4, 8).unwrap(),
            double.eval_with_normal(x, y, ny),
        ),
        (
            FmmC64::new_adjoint_double(adjoint, points.clone(), normals.clone(), 4, 8).unwrap(),
            adjoint.eval_with_normal(x, y, nx),
        ),
        (
            FmmC64::new_hypersingular(mixed, points, normals, 4, 8).unwrap(),
            mixed.eval_with_normals(x, y, nx, ny),
        ),
    ];
    for (engine, kernel) in engines {
        close(engine.direct(&q).unwrap()[0], kernel * q[1], 1.0e-14);
        close(engine.potentials(&q).unwrap()[0], kernel * q[1], 1.0e-14);
    }
}

#[test]
fn complex_storage_and_direct_envelopes_apply_before_pass_allocation_g0() {
    let points = (0..1024)
        .map(|j| {
            [
                (j % 16) as f64 / 16.0,
                ((j / 16) % 8) as f64 / 8.0,
                (j / 128) as f64 / 8.0,
            ]
        })
        .collect();
    let refused = FmmC64::new_hypersingular(
        HelmholtzHypersingular { k: 0.0 },
        points,
        vec![[1.0, 0.0, 0.0]; 1024],
        12,
        1,
    );
    assert!(matches!(
        refused,
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "complex coefficient-slot (two f64 per C64)",
            ..
        })
    ));
    let count = MAX_DIRECT_POINTS + 1;
    let points = (0..count)
        .map(|j| if j % 2 == 0 { [0.0; 3] } else { [1.0; 3] })
        .collect();
    let kernel = HelmholtzSingle { k: 0.0 };
    let fmm = FmmC64::new(&kernel, points, 2, 1000).unwrap();
    assert!(matches!(
        fmm.direct(&vec![C64::ONE; count]),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "complex direct pair",
            ..
        })
    ));
    assert_eq!(core::mem::size_of::<C64>(), 2 * core::mem::size_of::<f64>());
}

#[test]
fn declared_workers_leave_potentials_bitwise_g0() {
    // Several tree levels (leaf cap 2), all five operator shapes, with and
    // without the near-field cache, against the untouched serial path.
    let points: Vec<[f64; 3]> = (0..240)
        .map(|i| {
            let j = 0.07 * det::cos(i as f64);
            [
                (i % 8) as f64 + j,
                ((i / 8) % 6) as f64 - j,
                (i / 48) as f64 + 0.5 * j,
            ]
        })
        .collect();
    let normals: Vec<[f64; 3]> = points
        .iter()
        .map(|p| {
            let v = [p[0] + 1.0, p[1] + 2.0, p[2] + 3.0];
            let n = det::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            [v[0] / n, v[1] / n, v[2] / n]
        })
        .collect();
    let q: Vec<C64> = (0..240)
        .map(|i| C64::new(det::sin(i as f64), 0.2 - i as f64 / 240.0))
        .collect();
    let bits = |p: Vec<C64>| -> Vec<(u64, u64)> {
        p.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
    };
    let single = HelmholtzSingle { k: 0.6 };
    let gradient = HelmholtzGradient::new(0.6, 1).unwrap();
    let build = |shape: usize| {
        let (p, n) = (points.clone(), normals.clone());
        match shape {
            0 => FmmC64::new(&single, p, 3, 2),
            1 => FmmC64::new(&gradient, p, 3, 2),
            2 => FmmC64::new_double(HelmholtzDouble { k: 0.6 }, p, n, 3, 2),
            3 => FmmC64::new_adjoint_double(HelmholtzAdjointDouble { k: 0.6 }, p, n, 3, 2),
            _ => FmmC64::new_hypersingular(HelmholtzHypersingular { k: 0.6 }, p, n, 3, 2),
        }
        .unwrap()
    };
    for shape in 0..5 {
        for cached in [false, true] {
            let mut engine = build(shape);
            if cached {
                engine.cache_near_field(usize::MAX).unwrap();
            }
            assert!(engine.geometry.cells.keys().map(|k| k.0).max() >= Some(3));
            let serial = bits(engine.potentials(&q).unwrap());
            for workers in [2, 3, 7, 300] {
                engine.set_workers(workers).unwrap();
                assert_eq!(
                    bits(engine.potentials(&q).unwrap()),
                    serial,
                    "shape {shape}, cached {cached}, workers {workers}"
                );
            }
        }
    }
    // Each extra worker's M2L row is charged to the slot envelope.
    let mut engine = build(4);
    let base = engine.coefficient_slots;
    let row = 3usize.pow(3) * 3 * 3;
    engine.set_workers(5).unwrap();
    assert_eq!(engine.coefficient_slots, base + 4 * row);
    engine.set_workers(1).unwrap();
    assert_eq!(engine.coefficient_slots, base);
    assert!(matches!(
        engine.set_workers(usize::MAX / 4),
        Err(FmmError::WorkEnvelopeExceeded { .. })
    ));
    assert_eq!((engine.coefficient_slots, engine.workers), (base, 1));
    // A non-finite M2L block reports the serial path's first error.
    struct Broken;
    impl ComplexKernel for Broken {
        fn eval(&self, _: [f64; 3], _: [f64; 3]) -> C64 {
            C64::new(f64::NAN, 0.0)
        }
    }
    let mut broken = FmmC64::new(&Broken, points.clone(), 3, 2).unwrap();
    let serial = broken.potentials(&q).unwrap_err();
    broken.set_workers(4).unwrap();
    assert_eq!(broken.potentials(&q).unwrap_err(), serial);
}
