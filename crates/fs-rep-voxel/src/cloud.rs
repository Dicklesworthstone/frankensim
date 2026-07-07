//! Point clouds with normals — THE FITTING TARGET for scan-to-Region
//! workflows. Grid-hash accelerated k-NN/radius queries (brute-force
//! verified), PCA normal estimation (smallest covariance eigenvector via
//! cyclic Jacobi), and BFS orientation propagation over the k-NN graph.

use crate::VoxelError;
use fs_math::det;
use std::collections::BTreeMap;

/// A point cloud with optional per-point unit normals.
#[derive(Debug, Clone, PartialEq)]
pub struct PointCloud {
    /// Positions.
    pub points: Vec<[f64; 3]>,
    /// Unit normals (filled by [`PointCloud::estimate_normals`]).
    pub normals: Option<Vec<[f64; 3]>>,
    cell: f64,
    hash: BTreeMap<[i64; 3], Vec<usize>>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    det::sqrt(dot(a, a))
}

impl PointCloud {
    /// Build a cloud with a spatial hash of cell size `cell`.
    ///
    /// # Errors
    /// [`VoxelError::Cloud`] on an empty cloud, non-finite points, or a
    /// non-positive cell size.
    pub fn new(points: Vec<[f64; 3]>, cell: f64) -> Result<Self, VoxelError> {
        if points.is_empty() {
            return Err(VoxelError::Cloud {
                what: "empty point cloud".to_string(),
            });
        }
        if !(cell.is_finite() && cell > 0.0) {
            return Err(VoxelError::Cloud {
                what: format!("hash cell size {cell} must be positive"),
            });
        }
        if points.iter().flatten().any(|v| !v.is_finite()) {
            return Err(VoxelError::Cloud {
                what: "non-finite point coordinates".to_string(),
            });
        }
        let mut hash: BTreeMap<[i64; 3], Vec<usize>> = BTreeMap::new();
        for (i, p) in points.iter().enumerate() {
            hash.entry(Self::key(*p, cell)).or_default().push(i);
        }
        Ok(PointCloud {
            points,
            normals: None,
            cell,
            hash,
        })
    }

    fn key(p: [f64; 3], cell: f64) -> [i64; 3] {
        core::array::from_fn(|k| {
            #[allow(clippy::cast_possible_truncation)]
            {
                (p[k] / cell).floor() as i64
            }
        })
    }

    /// All indices within `radius` of `q`, sorted by (distance, index) —
    /// deterministic and brute-force verified.
    #[must_use]
    pub fn radius_query(&self, q: [f64; 3], radius: f64) -> Vec<usize> {
        let r2 = radius * radius;
        #[allow(clippy::cast_possible_truncation)]
        let reach = (radius / self.cell).ceil() as i64;
        let center = Self::key(q, self.cell);
        let mut hits: Vec<(f64, usize)> = Vec::new();
        for dx in -reach..=reach {
            for dy in -reach..=reach {
                for dz in -reach..=reach {
                    let cell_key = [center[0] + dx, center[1] + dy, center[2] + dz];
                    if let Some(bucket) = self.hash.get(&cell_key) {
                        for &i in bucket {
                            let d2 = {
                                let d = sub(self.points[i], q);
                                dot(d, d)
                            };
                            if d2 <= r2 {
                                hits.push((d2, i));
                            }
                        }
                    }
                }
            }
        }
        hits.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        hits.into_iter().map(|(_, i)| i).collect()
    }

    /// The `k` nearest neighbors of `q` (excluding exact self-matches is
    /// the CALLER's concern — indices sorted by distance then index).
    #[must_use]
    pub fn knn(&self, q: [f64; 3], k: usize) -> Vec<usize> {
        // Expand rings until enough candidates, then verify with one more
        // ring (grid-hash correctness: nearest may sit in the next ring).
        let mut radius = self.cell;
        loop {
            let hits = self.radius_query(q, radius);
            if hits.len() >= k || radius > self.cell * 1e6 {
                let confirm = self.radius_query(q, radius + self.cell);
                return confirm.into_iter().take(k).collect();
            }
            radius *= 2.0;
        }
    }

    /// Estimate unit normals from the `k`-neighborhood covariance
    /// (smallest eigenvector), then propagate a consistent orientation by
    /// BFS over the k-NN graph from the point with the largest z
    /// (its normal is seeded to point +z-ward).
    ///
    /// # Errors
    /// [`VoxelError::Cloud`] when `k < 3` or the cloud is smaller than
    /// `k + 1`.
    pub fn estimate_normals(&mut self, k: usize) -> Result<(), VoxelError> {
        if k < 3 || self.points.len() <= k {
            return Err(VoxelError::Cloud {
                what: format!("need k >= 3 and more than k={k} points"),
            });
        }
        let n = self.points.len();
        let mut normals = vec![[0.0f64; 3]; n];
        let mut neighborhoods = Vec::with_capacity(n);
        for i in 0..n {
            let nb: Vec<usize> = self
                .knn(self.points[i], k + 1)
                .into_iter()
                .filter(|&j| j != i)
                .take(k)
                .collect();
            let mut mean = [0.0f64; 3];
            for &j in &nb {
                for a in 0..3 {
                    mean[a] += self.points[j][a];
                }
            }
            #[allow(clippy::cast_precision_loss)]
            let inv = 1.0 / nb.len() as f64;
            for m in &mut mean {
                *m *= inv;
            }
            let mut cov = [[0.0f64; 3]; 3];
            for &j in &nb {
                let d = sub(self.points[j], mean);
                for a in 0..3 {
                    for b in 0..3 {
                        cov[a][b] += d[a] * d[b];
                    }
                }
            }
            normals[i] = smallest_eigenvector(cov);
            neighborhoods.push(nb);
        }
        // Orientation propagation: BFS from the topmost point, aligning
        // each normal with its parent's.
        let seed = (0..n)
            .max_by(|&a, &b| self.points[a][2].total_cmp(&self.points[b][2]))
            .expect("nonempty");
        if normals[seed][2] < 0.0 {
            normals[seed] = normals[seed].map(|v| -v);
        }
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::from([seed]);
        visited[seed] = true;
        loop {
            while let Some(i) = queue.pop_front() {
                for &j in &neighborhoods[i] {
                    if !visited[j] {
                        visited[j] = true;
                        if dot(normals[j], normals[i]) < 0.0 {
                            normals[j] = normals[j].map(|v| -v);
                        }
                        queue.push_back(j);
                    }
                }
            }
            // The k-NN graph can be DISCONNECTED (scan lines: dense along
            // a ring, sparse across rings). Restart each component from
            // its lowest-index point, aligning with the spatially nearest
            // already-oriented point so global consistency survives.
            let Some(restart) = (0..n).find(|&i| !visited[i]) else {
                break;
            };
            let mut radius = self.cell;
            let mut anchor = None;
            while anchor.is_none() && radius < self.cell * 1e9 {
                anchor = self
                    .radius_query(self.points[restart], radius)
                    .into_iter()
                    .find(|&j| visited[j]);
                radius *= 2.0;
            }
            match anchor {
                Some(j) => {
                    if dot(normals[restart], normals[j]) < 0.0 {
                        normals[restart] = normals[restart].map(|v| -v);
                    }
                }
                None => {
                    if normals[restart][2] < 0.0 {
                        normals[restart] = normals[restart].map(|v| -v);
                    }
                }
            }
            visited[restart] = true;
            queue.push_back(restart);
        }
        self.normals = Some(normals);
        Ok(())
    }
}

/// Unit eigenvector of the smallest eigenvalue of a symmetric 3×3
/// matrix, via cyclic Jacobi (deterministic sweep order).
fn smallest_eigenvector(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        let off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if off < 1e-30 {
            break;
        }
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            if a[p][q].abs() < 1e-300 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
            let t = theta.signum() / (theta.abs() + det::sqrt(theta * theta + 1.0));
            let c = 1.0 / det::sqrt(t * t + 1.0);
            let s = t * c;
            for k in 0..3 {
                let (apk, aqk) = (a[p][k], a[q][k]);
                a[p][k] = c * apk - s * aqk;
                a[q][k] = s * apk + c * aqk;
            }
            for k in 0..3 {
                let (akp, akq) = (a[k][p], a[k][q]);
                a[k][p] = c * akp - s * akq;
                a[k][q] = s * akp + c * akq;
                let (vkp, vkq) = (v[k][p], v[k][q]);
                v[k][p] = c * vkp - s * vkq;
                v[k][q] = s * vkp + c * vkq;
            }
        }
    }
    let evs = [a[0][0], a[1][1], a[2][2]];
    let smallest = (0..3)
        .min_by(|&i, &j| evs[i].total_cmp(&evs[j]))
        .expect("three eigenvalues");
    let col = [v[0][smallest], v[1][smallest], v[2][smallest]];
    let len = norm(col).max(f64::MIN_POSITIVE);
    col.map(|x| x / len)
}
