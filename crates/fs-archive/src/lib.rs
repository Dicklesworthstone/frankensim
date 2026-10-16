//! fs-archive — quality-diversity archives (MAP-Elites / CVT). Layer: L4.
//!
//! Optimization returns ONE point. ILLUMINATION returns a map: the best solution
//! found in every behavioral niche. That is the quality-diversity idea, and it
//! is exactly what a design study wants — not "the single stiffest bracket" but
//! "the best bracket at each mass, each footprint, each lip topology."
//!
//! [`MapElites`] discretizes a behavior-descriptor space into a grid and keeps
//! one ELITE per cell (the highest-fitness solution with that descriptor); a
//! worse candidate in a filled cell is rejected, a better one replaces it. So
//! [`MapElites::coverage`] (fraction of niches filled) and
//! [`MapElites::qd_score`] (total elite fitness) only ever rise. [`CvtArchive`]
//! is the same illumination over a Centroidal-Voronoi tessellation (assign to
//! the nearest centroid) for high-dimensional descriptors, and [`novelty`]
//! scores how far a descriptor is from the current archive (exploration
//! pressure). Deterministic; no dependencies.

use std::collections::BTreeMap;

/// An elite: the best solution recorded for its behavioral niche.
#[derive(Debug, Clone, PartialEq)]
pub struct Elite {
    /// The solution (genotype / design vector).
    pub solution: Vec<f64>,
    /// Its behavior descriptor.
    pub descriptor: Vec<f64>,
    /// Its fitness (higher is better).
    pub fitness: f64,
}

/// A MAP-Elites illumination archive over a gridded behavior space.
#[derive(Debug, Clone)]
pub struct MapElites {
    lo: Vec<f64>,
    hi: Vec<f64>,
    bins: Vec<usize>,
    cells: BTreeMap<Vec<usize>, Elite>,
}

impl MapElites {
    /// A new archive over `[lo, hi]` per descriptor dimension, discretized into
    /// `bins` cells per dimension.
    ///
    /// # Panics
    /// If dimensions disagree, any bin count is zero, bounds are non-finite,
    /// or any upper bound is not greater than its lower bound.
    #[must_use]
    pub fn new(lo: Vec<f64>, hi: Vec<f64>, bins: Vec<usize>) -> MapElites {
        assert!(
            lo.len() == hi.len() && lo.len() == bins.len(),
            "descriptor dims disagree"
        );
        assert!(
            !lo.is_empty() && bins.iter().all(|&b| b > 0),
            "need >=1 dim and positive bins"
        );
        assert!(
            lo.iter().zip(&hi).all(|(l, h)| h > l),
            "each hi must exceed lo"
        );
        assert!(
            lo.iter().chain(&hi).all(|v| v.is_finite()),
            "bounds must be finite"
        );
        MapElites {
            lo,
            hi,
            bins,
            cells: BTreeMap::new(),
        }
    }

    /// The discrete cell index of a descriptor (clamped into the grid).
    ///
    /// # Panics
    /// If the descriptor dimension disagrees with the archive or contains a
    /// non-finite value.
    #[must_use]
    pub fn cell_of(&self, descriptor: &[f64]) -> Vec<usize> {
        assert_descriptor("descriptor", descriptor, self.lo.len());
        descriptor
            .iter()
            .zip(&self.lo)
            .zip(&self.hi)
            .zip(&self.bins)
            .map(|(((&d, &l), &h), &b)| {
                let frac = (d - l) / (h - l);
                let idx = (frac * b as f64).floor();
                // clamp into [0, b-1].
                if idx < 0.0 {
                    0
                } else if idx >= b as f64 {
                    b - 1
                } else {
                    idx as usize
                }
            })
            .collect()
    }

    /// Try to add a solution. Returns `true` if it became an elite (a new niche
    /// or a strict improvement over the incumbent).
    ///
    /// # Panics
    /// If the descriptor dimension disagrees with the archive or fitness is
    /// negative / non-finite.
    pub fn add(&mut self, solution: Vec<f64>, descriptor: Vec<f64>, fitness: f64) -> bool {
        assert_non_negative_fitness(fitness);
        let cell = self.cell_of(&descriptor);
        let improve = self.cells.get(&cell).is_none_or(|e| fitness > e.fitness);
        if improve {
            self.cells.insert(
                cell,
                Elite {
                    solution,
                    descriptor,
                    fitness,
                },
            );
        }
        improve
    }

    /// The total number of cells in the grid.
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.bins.iter().product()
    }

    /// The number of filled niches.
    #[must_use]
    pub fn num_elites(&self) -> usize {
        self.cells.len()
    }

    /// The fraction of niches filled.
    #[must_use]
    pub fn coverage(&self) -> f64 {
        self.num_elites() as f64 / self.capacity() as f64
    }

    /// The quality-diversity score: the total fitness over all elites.
    #[must_use]
    pub fn qd_score(&self) -> f64 {
        self.cells.values().map(|e| e.fitness).sum()
    }

    /// The elite in the cell of a descriptor, if any.
    ///
    /// # Panics
    /// If the descriptor dimension disagrees with the archive or contains a
    /// non-finite value.
    #[must_use]
    pub fn elite_at(&self, descriptor: &[f64]) -> Option<&Elite> {
        self.cells.get(&self.cell_of(descriptor))
    }

    /// The single highest-fitness elite.
    #[must_use]
    pub fn best(&self) -> Option<&Elite> {
        self.cells.values().max_by(|a, b| {
            a.fitness
                .partial_cmp(&b.fitness)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// All elites (deterministic cell order).
    pub fn elites(&self) -> impl Iterator<Item = &Elite> {
        self.cells.values()
    }
}

/// A CVT (Centroidal Voronoi Tessellation) archive: illumination over
/// arbitrary centroids, assigning each descriptor to its nearest centroid.
#[derive(Debug, Clone)]
pub struct CvtArchive {
    centroids: Vec<Vec<f64>>,
    dim: usize,
    cells: BTreeMap<usize, Elite>,
}

impl CvtArchive {
    /// A new archive over the given centroids.
    ///
    /// # Panics
    /// If `centroids` is empty, zero-dimensional, non-finite, or not all the
    /// same dimension.
    #[must_use]
    pub fn new(centroids: Vec<Vec<f64>>) -> CvtArchive {
        assert!(!centroids.is_empty(), "need at least one centroid");
        let dim = centroids.first().map_or(0, Vec::len);
        assert!(dim > 0, "centroids need at least one dimension");
        for centroid in &centroids {
            assert_descriptor("centroid", centroid, dim);
        }
        CvtArchive {
            centroids,
            dim,
            cells: BTreeMap::new(),
        }
    }

    /// The index of the nearest centroid to a descriptor.
    ///
    /// # Panics
    /// If the descriptor dimension disagrees with the centroids or contains a
    /// non-finite value.
    #[must_use]
    pub fn nearest_centroid(&self, descriptor: &[f64]) -> usize {
        assert_descriptor("descriptor", descriptor, self.dim);
        self.centroids
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| compare_distances(a, b, descriptor))
            .map_or(0, |(i, _)| i)
    }

    /// Try to add a solution (nearest-centroid niche). Returns `true` if it
    /// became an elite.
    ///
    /// # Panics
    /// If the descriptor dimension disagrees with the centroids or fitness is
    /// negative / non-finite.
    pub fn add(&mut self, solution: Vec<f64>, descriptor: Vec<f64>, fitness: f64) -> bool {
        assert_non_negative_fitness(fitness);
        let cell = self.nearest_centroid(&descriptor);
        let improve = self.cells.get(&cell).is_none_or(|e| fitness > e.fitness);
        if improve {
            self.cells.insert(
                cell,
                Elite {
                    solution,
                    descriptor,
                    fitness,
                },
            );
        }
        improve
    }

    /// The number of centroids (niche capacity).
    #[must_use]
    pub fn capacity(&self) -> usize {
        self.centroids.len()
    }

    /// The number of filled niches.
    #[must_use]
    pub fn num_elites(&self) -> usize {
        self.cells.len()
    }

    /// The fraction of niches filled.
    #[must_use]
    pub fn coverage(&self) -> f64 {
        self.num_elites() as f64 / self.capacity() as f64
    }

    /// The quality-diversity score (total elite fitness).
    #[must_use]
    pub fn qd_score(&self) -> f64 {
        self.cells.values().map(|e| e.fitness).sum()
    }

    /// The single highest-fitness elite.
    #[must_use]
    pub fn best(&self) -> Option<&Elite> {
        self.cells.values().max_by(|a, b| {
            a.fitness
                .partial_cmp(&b.fitness)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
}

/// The NOVELTY of a descriptor relative to a set of others: the mean Euclidean
/// distance to its `k` nearest neighbours (exploration pressure). An empty
/// neighbour set or a request for zero neighbours is maximally novel (`+∞`).
/// A nonempty neighbour set is validated even when `k == 0`.
///
/// # Panics
/// If `others` is non-empty and descriptor dimensions disagree, are zero, or
/// contain non-finite values.
#[must_use]
pub fn novelty(descriptor: &[f64], others: &[Vec<f64>], k: usize) -> f64 {
    if others.is_empty() {
        return f64::INFINITY;
    }
    assert_descriptor("descriptor", descriptor, descriptor.len());
    for other in others {
        assert_descriptor("neighbour descriptor", other, descriptor.len());
    }
    if k == 0 {
        return f64::INFINITY;
    }
    let mut dists: Vec<DistanceMagnitude> = others
        .iter()
        .map(|other| distance_magnitude(other, descriptor))
        .collect();
    dists.sort_by(|a, b| a.compare(*b));
    let kk = k.min(dists.len());
    mean_distances(&dists[..kk])
}

fn assert_descriptor(label: &str, values: &[f64], expected_dim: usize) {
    assert!(expected_dim > 0, "{label} needs at least one dimension");
    assert!(
        values.len() == expected_dim,
        "{label} dimension {} does not match expected dimension {expected_dim}",
        values.len()
    );
    assert!(
        values.iter().all(|v| v.is_finite()),
        "{label} must be finite"
    );
}

fn assert_non_negative_fitness(fitness: f64) {
    assert!(
        fitness.is_finite() && fitness >= 0.0,
        "fitness must be finite and non-negative"
    );
}

fn squared_distance(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

#[derive(Debug, Clone, Copy)]
struct DistanceMagnitude {
    squared: f64,
    scale: f64,
    factor: f64,
}

impl DistanceMagnitude {
    const ZERO: DistanceMagnitude = DistanceMagnitude {
        squared: 0.0,
        scale: 0.0,
        factor: 0.0,
    };

    fn materialize(self) -> f64 {
        self.scale * self.factor
    }

    fn compare(self, other: DistanceMagnitude) -> std::cmp::Ordering {
        match (self.squared.is_finite(), other.squared.is_finite()) {
            (true, true) => self
                .squared
                .partial_cmp(&other.squared)
                .unwrap_or(std::cmp::Ordering::Equal),
            _ => {
                if self.scale >= other.scale {
                    self.factor
                        .partial_cmp(&(other.factor * (other.scale / self.scale)))
                        .unwrap_or(std::cmp::Ordering::Equal)
                } else {
                    (self.factor * (self.scale / other.scale))
                        .partial_cmp(&other.factor)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }
            }
        }
    }

    fn ratio_to(self, maximum: DistanceMagnitude) -> f64 {
        if self.scale == 0.0 {
            return 0.0;
        }
        let ratio = if self.scale <= maximum.scale {
            (self.factor * (self.scale / maximum.scale)) / maximum.factor
        } else {
            self.factor / (maximum.factor * (maximum.scale / self.scale))
        };
        ratio.min(1.0)
    }
}

fn distance_magnitude(a: &[f64], b: &[f64]) -> DistanceMagnitude {
    let squared = squared_distance(a, b);
    if squared.is_finite() {
        let distance = squared.sqrt();
        return if distance == 0.0 {
            DistanceMagnitude::ZERO
        } else {
            DistanceMagnitude {
                squared,
                scale: distance,
                factor: 1.0,
            }
        };
    }
    let scale = a
        .iter()
        .zip(b)
        .map(|(x, y)| x.abs().max(y.abs()))
        .fold(0.0_f64, f64::max);
    let normalized_squared: f64 = a
        .iter()
        .zip(b)
        .map(|(x, y)| {
            let delta = x / scale - y / scale;
            delta * delta
        })
        .sum();
    DistanceMagnitude {
        squared,
        scale,
        factor: normalized_squared.sqrt(),
    }
}

fn compare_distances(a: &[f64], b: &[f64], descriptor: &[f64]) -> std::cmp::Ordering {
    let a_squared = squared_distance(a, descriptor);
    let b_squared = squared_distance(b, descriptor);
    match (a_squared.is_finite(), b_squared.is_finite()) {
        (true, true) => a_squared
            .partial_cmp(&b_squared)
            .unwrap_or(std::cmp::Ordering::Equal),
        _ => distance_magnitude(a, descriptor).compare(distance_magnitude(b, descriptor)),
    }
}

fn mean_distances(values: &[DistanceMagnitude]) -> f64 {
    debug_assert!(!values.is_empty());
    let count = values.len() as f64;
    let direct_sum = values.iter().map(|value| value.materialize()).sum::<f64>();
    if direct_sum.is_finite() {
        return direct_sum / count;
    }
    let maximum = values
        .iter()
        .copied()
        .max_by(|a, b| a.compare(*b))
        .expect("nonempty distance slice");
    if maximum.scale == 0.0 {
        return 0.0;
    }
    let ratio_sum = values
        .iter()
        .map(|value| value.ratio_to(maximum))
        .sum::<f64>();
    maximum.scale * (maximum.factor * (ratio_sum / count))
}
