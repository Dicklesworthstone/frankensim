//! Battery for quality-diversity archives (fs-archive). Covers descriptor
//! discretization, the MAP-Elites keep-the-best-per-niche rule, the
//! quality-diversity property (diversity preserved, not just the global best),
//! coverage + QD-score, the CVT archive, novelty scoring, and determinism.

use fs_archive::{CvtArchive, MapElites, novelty};

fn grid() -> MapElites {
    MapElites::new(vec![0.0, 0.0], vec![1.0, 1.0], vec![4, 4])
}

fn panics(f: impl FnOnce() + std::panic::UnwindSafe) -> bool {
    std::panic::catch_unwind(f).is_err()
}

#[test]
fn descriptors_map_to_the_right_cells() {
    let a = grid();
    assert_eq!(a.cell_of(&[0.1, 0.1]), vec![0, 0]);
    assert_eq!(a.cell_of(&[0.9, 0.9]), vec![3, 3]);
    assert_eq!(a.cell_of(&[0.5, 0.5]), vec![2, 2]);
    // out-of-range descriptors clamp into the grid.
    assert_eq!(a.cell_of(&[1.5, -0.5]), vec![3, 0]);
}

#[test]
fn map_elites_keeps_the_best_per_niche() {
    let mut a = grid();
    assert!(a.add(vec![1.0], vec![0.2, 0.2], 5.0)); // new niche -> inserted
    assert_eq!(a.num_elites(), 1);
    // a worse candidate in the same cell is rejected.
    assert!(!a.add(vec![2.0], vec![0.21, 0.21], 3.0));
    assert!((a.elite_at(&[0.2, 0.2]).unwrap().fitness - 5.0).abs() < 1e-12);
    // a better candidate replaces the incumbent.
    assert!(a.add(vec![3.0], vec![0.22, 0.22], 9.0));
    assert!((a.elite_at(&[0.2, 0.2]).unwrap().fitness - 9.0).abs() < 1e-12);
    assert_eq!(a.num_elites(), 1); // still one niche
}

#[test]
fn illumination_preserves_diversity_not_just_the_best() {
    let mut a = grid();
    a.add(vec![1.0], vec![0.1, 0.1], 10.0); // high fitness, niche X
    a.add(vec![2.0], vec![0.9, 0.9], 1.0); // LOW fitness, niche Y
    // both are kept -> a pure optimizer would have discarded the low-fitness one.
    assert_eq!(a.num_elites(), 2);
    assert!((a.best().unwrap().fitness - 10.0).abs() < 1e-12);
    assert!((a.qd_score() - 11.0).abs() < 1e-12);
    assert!((a.coverage() - 2.0 / 16.0).abs() < 1e-12);
}

#[test]
fn coverage_and_qd_score_are_monotone() {
    let mut a = grid();
    let (mut cov, mut qd) = (0.0, 0.0);
    for i in 0..12 {
        let x = f64::from(i) / 12.0;
        a.add(vec![x], vec![x, 1.0 - x], f64::from(i));
        assert!(a.coverage() >= cov && a.qd_score() >= qd);
        cov = a.coverage();
        qd = a.qd_score();
    }
    assert!(a.num_elites() > 1);
}

#[test]
fn the_cvt_archive_assigns_to_the_nearest_centroid() {
    let mut c = CvtArchive::new(vec![vec![0.0, 0.0], vec![1.0, 1.0], vec![0.0, 1.0]]);
    assert_eq!(c.nearest_centroid(&[0.1, 0.1]), 0);
    assert_eq!(c.nearest_centroid(&[0.9, 0.8]), 1);
    assert_eq!(c.nearest_centroid(&[0.1, 0.9]), 2);
    c.add(vec![1.0], vec![0.1, 0.1], 4.0);
    c.add(vec![2.0], vec![0.9, 0.9], 7.0);
    assert!(!c.add(vec![3.0], vec![0.15, 0.05], 2.0)); // worse in centroid 0
    assert_eq!(c.num_elites(), 2);
    assert!((c.qd_score() - 11.0).abs() < 1e-12);
    assert!((c.coverage() - 2.0 / 3.0).abs() < 1e-12);
}

#[test]
fn novelty_rewards_distance_from_the_archive() {
    let others = vec![vec![0.0, 0.0], vec![0.0, 1.0], vec![1.0, 0.0]];
    let far = novelty(&[5.0, 5.0], &others, 2);
    let near = novelty(&[0.0, 0.5], &others, 2);
    assert!(far > near, "far {far} should exceed near {near}");
    // an empty archive is maximally novel.
    assert!(novelty(&[0.0], &[], 3).is_infinite());
    // Asking for no neighbours is also maximally novel after validation.
    assert!(novelty(&[0.0], &[vec![0.0]], 0).is_infinite());
}

#[test]
fn malformed_dimensions_and_fitness_are_rejected() {
    assert!(panics(|| {
        let _ = grid().cell_of(&[0.2]);
    }));
    assert!(panics(|| {
        let mut a = grid();
        let _ = a.add(vec![1.0], vec![0.2, 0.2], -1.0);
    }));
    assert!(panics(|| {
        let _ = CvtArchive::new(vec![vec![0.0, 0.0], vec![1.0]]);
    }));
    assert!(panics(|| {
        let c = CvtArchive::new(vec![vec![0.0, 0.0]]);
        let _ = c.nearest_centroid(&[0.0]);
    }));
    assert!(panics(|| {
        let _ = novelty(&[0.0, 0.0], &[vec![0.0]], 1);
    }));
    // k=0 is an output policy, not a validation bypass for a nonempty archive.
    assert!(panics(|| {
        let _ = novelty(&[], &[vec![0.0]], 0);
    }));
    assert!(panics(|| {
        let _ = novelty(&[f64::NAN], &[vec![0.0]], 0);
    }));
    assert!(panics(|| {
        let _ = novelty(&[0.0], &[vec![f64::INFINITY]], 0);
    }));
    assert!(panics(|| {
        let _ = novelty(&[0.0, 0.0], &[vec![0.0]], 0);
    }));
}

#[test]
fn archives_are_deterministic() {
    let build = || {
        let mut a = grid();
        for i in 0..20 {
            let x = f64::from(i % 4) / 4.0 + 0.05;
            let y = f64::from(i % 3) / 3.0 + 0.05;
            a.add(vec![f64::from(i)], vec![x, y], f64::from(i));
        }
        a
    };
    let (a, b) = (build(), build());
    assert_eq!(a.num_elites(), b.num_elites());
    assert_eq!(a.qd_score().to_bits(), b.qd_score().to_bits());
    assert_eq!(
        a.best().unwrap().fitness.to_bits(),
        b.best().unwrap().fitness.to_bits()
    );
}
