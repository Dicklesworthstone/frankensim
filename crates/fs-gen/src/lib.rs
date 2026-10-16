//! fs-gen — PROPOSAL-ONLY generative models (patch Rev N, bead 7tv.19;
//! [M] — behind the `proposal-gen` feature). The rule is absolute:
//! generative models may PROPOSE (initial designs, restart seeds,
//! mutation directions); they may NEVER CERTIFY (physics outcomes,
//! safety, feasibility, final rankings). The COMPILER enforces the
//! epistemics: every generator output is a [`Proposal<T>`] whose
//! payload is PRIVATE — the only exit is [`Proposal::promote`], which
//! runs the caller's validation machinery. There is no Evidence field
//! on a proposal because a proposal IS NOT EVIDENCE.
//!
//! v0 generators are honest classical density machinery fitted on the
//! ledger corpus (the system learns its OWN design distribution);
//! diffusion/flow and transformer program proposers are the documented
//! [M] growth path, not smuggled in.
#![cfg(feature = "proposal-gen")]

/// The provenance card every proposal carries (training-data
/// provenance from ledger hashes; determinism class documented).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelCard {
    /// Generative-model identity/version.
    pub model: String,
    /// Content hash of the training corpus (ledger provenance).
    pub corpus_hash: String,
    /// Determinism class ("counter-based: bit-replayable per seed").
    pub determinism: String,
}

/// A PROPOSAL: a payload the type system quarantines. No Evidence
/// fields, no public payload access — the ONLY way out is
/// [`Proposal::promote`] through a validator.
///
/// ```compile_fail
/// // The payload is private: proposals cannot leak into certified
/// // paths without promotion. THIS MUST NOT COMPILE:
/// let p = fs_gen::Proposal::new(vec![1.0], fs_gen::ModelCard {
///     model: "m".into(), corpus_hash: "h".into(), determinism: "d".into(),
/// });
/// let sneak: &Vec<f64> = &p.payload;
/// ```
#[derive(Debug, Clone)]
pub struct Proposal<T> {
    payload: T,
    /// The provenance card (public — provenance is not evidence).
    pub card: ModelCard,
}

/// A validator's structured rejection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rejected {
    /// Why the proposal failed validation.
    pub reason: String,
}

impl<T> Proposal<T> {
    /// Wrap a generator output (generators only; downstream code never
    /// constructs proposals).
    #[must_use]
    pub fn new(payload: T, card: ModelCard) -> Proposal<T> {
        Proposal { payload, card }
    }

    /// THE ONLY EXIT: run the standard validation machinery. On
    /// success the payload leaves quarantine (the caller then attaches
    /// whatever Evidence its validation actually earned); on failure a
    /// structured rejection comes back and the payload stays inside.
    ///
    /// # Errors
    /// [`Rejected`] with the validator's reason.
    pub fn promote<E>(self, validate: impl FnOnce(&T) -> Result<(), E>) -> Result<T, Rejected>
    where
        E: std::fmt::Display,
    {
        match validate(&self.payload) {
            Ok(()) => Ok(self.payload),
            Err(e) => Err(Rejected {
                reason: e.to_string(),
            }),
        }
    }

    /// Peek for LOGGING ONLY: a read-only view gated behind an
    /// explicitly named method so audits can see what was proposed —
    /// the name is the warning.
    #[must_use]
    pub fn inspect_for_logging_only(&self) -> &T {
        &self.payload
    }
}

/// Counter-based uniform in [0,1) (splitmix-style; bit-replayable).
fn unit(seed: u64, k: u64) -> f64 {
    let mut z = seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k.wrapping_add(1));
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// Approximate standard normal (12-uniform sum; deterministic).
fn gauss(seed: u64, k: u64) -> f64 {
    (0..12)
        .map(|j| unit(seed, k.wrapping_mul(12).wrapping_add(j)))
        .sum::<f64>()
        - 6.0
}

/// FNV-1a corpus hash (ledger-style content address).
fn corpus_hash(corpus: &[Vec<f64>]) -> String {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for b in (corpus.len() as u64).to_le_bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    for row in corpus {
        for b in (row.len() as u64).to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        for v in row {
            for b in v.to_bits().to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
    }
    format!("{h:016x}")
}

fn corpus_dimension(corpus: &[Vec<f64>], label: &str) -> usize {
    assert!(!corpus.is_empty(), "{label} needs a corpus");
    let d = corpus[0].len();
    assert!(d > 0, "{label} needs non-empty design vectors");
    assert!(
        corpus
            .iter()
            .all(|row| row.len() == d && row.iter().all(|v| v.is_finite())),
        "{label} corpus rows must share one finite, non-empty dimension"
    );
    d
}

/// A CORPUS-FITTED SHAPE PRIOR: kernel density estimate over the
/// ledger's design vectors — the simplest honest "the system learns
/// its own design distribution". Sampling = pick a corpus anchor,
/// jitter by the fitted bandwidth (counter-based, replayable).
#[derive(Debug, Clone)]
pub struct ShapePrior {
    anchors: Vec<Vec<f64>>,
    bandwidth: f64,
    card: ModelCard,
}

impl ShapePrior {
    /// Fit on the corpus with Silverman-style bandwidth.
    ///
    /// # Panics
    /// On an empty corpus, zero-dimensional rows, non-finite values, or
    /// inconsistent row dimensions.
    #[must_use]
    pub fn fit(corpus: &[Vec<f64>]) -> ShapePrior {
        let d = corpus_dimension(corpus, "shape prior");
        // Mean per-coordinate std for the bandwidth scale.
        let mut var_sum = 0.0f64;
        for j in 0..d {
            #[allow(clippy::cast_precision_loss)]
            let mean = corpus.iter().map(|r| r[j]).sum::<f64>() / corpus.len() as f64;
            #[allow(clippy::cast_precision_loss)]
            let var =
                corpus.iter().map(|r| (r[j] - mean).powi(2)).sum::<f64>() / corpus.len() as f64;
            var_sum += var;
        }
        #[allow(clippy::cast_precision_loss)]
        let sigma = (var_sum / d as f64).sqrt().max(1e-9);
        #[allow(clippy::cast_precision_loss)]
        let bandwidth = sigma * (corpus.len() as f64).powf(-1.0 / (d as f64 + 4.0));
        ShapePrior {
            anchors: corpus.to_vec(),
            bandwidth,
            card: ModelCard {
                model: "kde-shape-prior-v0".to_string(),
                corpus_hash: corpus_hash(corpus),
                determinism: "counter-based: bit-replayable per seed".to_string(),
            },
        }
    }

    /// The prior density (up to normalization) at a point — the
    /// "prior is confident" signal acquisition uses.
    #[must_use]
    pub fn density(&self, x: &[f64]) -> f64 {
        assert_eq!(
            x.len(),
            self.anchors[0].len(),
            "density query dimension must match the fitted corpus"
        );
        assert!(
            x.iter().all(|v| v.is_finite()),
            "density query must be finite"
        );
        let h2 = 2.0 * self.bandwidth * self.bandwidth;
        self.anchors
            .iter()
            .map(|a| {
                let d2: f64 = a.iter().zip(x).map(|(ai, xi)| (ai - xi).powi(2)).sum();
                (-d2 / h2).exp()
            })
            .sum()
    }

    /// Sample one PROPOSAL (anchor + bandwidth jitter).
    #[must_use]
    pub fn propose(&self, seed: u64) -> Proposal<Vec<f64>> {
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let idx = (unit(seed, 0) * self.anchors.len() as f64) as usize;
        let anchor = &self.anchors[idx.min(self.anchors.len() - 1)];
        let sample: Vec<f64> = anchor
            .iter()
            .enumerate()
            .map(|(j, a)| a + self.bandwidth * gauss(seed, 1 + j as u64))
            .collect();
        Proposal::new(sample, self.card.clone())
    }
}

/// A COVARIANCE-SHAPED MUTATION KERNEL: perturbation directions shaped
/// by the corpus covariance's dominant axes (power iteration) — the
/// learned alternative to isotropic noise for level-set mutations.
#[derive(Debug, Clone)]
pub struct MutationKernel {
    /// Dominant corpus direction (unit).
    principal: Vec<f64>,
    /// Std along the principal direction vs isotropic floor.
    along: f64,
    /// Isotropic floor.
    floor: f64,
    card: ModelCard,
}

impl MutationKernel {
    /// Fit from the corpus (dominant covariance axis by power
    /// iteration — deterministic start).
    ///
    /// # Panics
    /// On an empty corpus, zero-dimensional rows, non-finite values, or
    /// inconsistent row dimensions.
    #[must_use]
    pub fn fit(corpus: &[Vec<f64>]) -> MutationKernel {
        let d = corpus_dimension(corpus, "mutation kernel");
        #[allow(clippy::cast_precision_loss)]
        let n = corpus.len() as f64;
        let mean: Vec<f64> = (0..d)
            .map(|j| corpus.iter().map(|r| r[j]).sum::<f64>() / n)
            .collect();
        // Power iteration on the covariance (matrix-free).
        let mut v: Vec<f64> = (0..d).map(|j| 1.0 + 0.1 * (j as f64)).collect();
        for _ in 0..60 {
            let mut next = vec![0.0f64; d];
            for row in corpus {
                let dot: f64 = row
                    .iter()
                    .zip(&mean)
                    .zip(&v)
                    .map(|((r, m), vi)| (r - m) * vi)
                    .sum();
                for (nj, (rj, mj)) in next.iter_mut().zip(row.iter().zip(&mean)) {
                    *nj += dot * (rj - mj) / n;
                }
            }
            let norm = next.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm <= 1e-24 {
                v.fill(0.0);
                v[0] = 1.0;
                break;
            }
            v = next.into_iter().map(|x| x / norm).collect();
        }
        // Variance along the axis + mean residual variance.
        let along_var: f64 = corpus
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&mean)
                    .zip(&v)
                    .map(|((r, m), vi)| (r - m) * vi)
                    .sum::<f64>()
                    .powi(2)
            })
            .sum::<f64>()
            / n;
        #[allow(clippy::cast_precision_loss)]
        let total_var: f64 = corpus
            .iter()
            .map(|row| {
                row.iter()
                    .zip(&mean)
                    .map(|(r, m)| (r - m).powi(2))
                    .sum::<f64>()
            })
            .sum::<f64>()
            / n;
        let residual = ((total_var - along_var) / (d as f64 - 1.0).max(1.0)).max(1e-12);
        MutationKernel {
            principal: v,
            along: along_var.sqrt(),
            floor: residual.sqrt(),
            card: ModelCard {
                model: "cov-mutation-kernel-v0".to_string(),
                corpus_hash: corpus_hash(corpus),
                determinism: "counter-based: bit-replayable per seed".to_string(),
            },
        }
    }

    /// Propose a mutation DIRECTION (scaled): principal-heavy jitter.
    #[must_use]
    pub fn propose(&self, seed: u64) -> Proposal<Vec<f64>> {
        let d = self.principal.len();
        let z_along = gauss(seed, 0) * self.along;
        let dir: Vec<f64> = (0..d)
            .map(|j| z_along * self.principal[j] + self.floor * gauss(seed, 1 + j as u64))
            .collect();
        Proposal::new(dir, self.card.clone())
    }
}

/// A TRUSS-GRAPH candidate generator: edges sampled from the corpus's
/// empirical degree bias (FrankenNetworkx-class priors beyond k-NN —
/// v0 keeps the graph as an edge list).
#[derive(Debug, Clone)]
pub struct GraphGenerator {
    n_nodes: usize,
    /// Empirical per-node attachment weight from the corpus.
    weights: Vec<f64>,
    card: ModelCard,
}

impl GraphGenerator {
    /// Fit attachment weights from corpus edge lists.
    #[must_use]
    pub fn fit(n_nodes: usize, corpus_edges: &[Vec<(usize, usize)>]) -> GraphGenerator {
        assert!(n_nodes > 0, "graph generator needs at least one node");
        let mut weights = vec![1.0f64; n_nodes];
        for graph in corpus_edges {
            for &(a, b) in graph {
                if a < n_nodes {
                    weights[a] += 1.0;
                }
                if b < n_nodes {
                    weights[b] += 1.0;
                }
            }
        }
        let flat: Vec<f64> = corpus_edges
            .iter()
            .flat_map(|g| g.iter().flat_map(|&(a, b)| [a as f64, b as f64]))
            .collect();
        GraphGenerator {
            n_nodes,
            weights,
            card: ModelCard {
                model: "degree-bias-graph-v0".to_string(),
                corpus_hash: corpus_hash(&[flat]),
                determinism: "counter-based: bit-replayable per seed".to_string(),
            },
        }
    }

    /// Propose a candidate edge list (`m` edges, degree-biased,
    /// self-loop-free, deterministic per seed).
    #[must_use]
    pub fn propose(&self, seed: u64, m: usize) -> Proposal<Vec<(usize, usize)>> {
        let max_edges = self.n_nodes.saturating_mul(self.n_nodes.saturating_sub(1)) / 2;
        assert!(
            m <= max_edges,
            "requested edge count exceeds the simple undirected graph capacity"
        );
        let total: f64 = self.weights.iter().sum();
        let pick = |u: f64| -> usize {
            let mut acc = 0.0;
            for (i, w) in self.weights.iter().enumerate() {
                acc += w / total;
                if u <= acc {
                    return i;
                }
            }
            self.n_nodes - 1
        };
        let mut edges = Vec::with_capacity(m);
        let mut k = 0u64;
        while edges.len() < m && k < 20 * m as u64 {
            let a = pick(unit(seed, 2 * k));
            let b = pick(unit(seed, 2 * k + 1));
            k += 1;
            if a != b && !edges.contains(&(a.min(b), a.max(b))) {
                edges.push((a.min(b), a.max(b)));
            }
        }
        for a in 0..self.n_nodes {
            for b in (a + 1)..self.n_nodes {
                if edges.len() == m {
                    break;
                }
                if !edges.contains(&(a, b)) {
                    edges.push((a, b));
                }
            }
        }
        Proposal::new(edges, self.card.clone())
    }
}

/// ACTIVE-LEARNING acquisition: propose where the ARCHIVE is sparse
/// AND the PRIOR is confident — exploration with taste. Returns the
/// top-k candidate proposals ranked by `density × min-dist-to-archive`.
#[must_use]
pub fn acquire(
    prior: &ShapePrior,
    archive: &[Vec<f64>],
    n_candidates: usize,
    k: usize,
    seed: u64,
) -> Vec<Proposal<Vec<f64>>> {
    let dim = prior.anchors[0].len();
    assert!(
        archive
            .iter()
            .all(|row| row.len() == dim && row.iter().all(|v| v.is_finite())),
        "archive rows must match the prior dimension and be finite"
    );
    let mut scored: Vec<(f64, Proposal<Vec<f64>>)> = (0..n_candidates)
        .map(|i| {
            let p = prior.propose(seed.wrapping_add((i as u64).wrapping_mul(7919)));
            let x = p.inspect_for_logging_only();
            let dens = prior.density(x);
            let dist = archive
                .iter()
                .map(|a| {
                    a.iter()
                        .zip(x)
                        .map(|(ai, xi)| (ai - xi).powi(2))
                        .sum::<f64>()
                        .sqrt()
                })
                .fold(f64::INFINITY, f64::min)
                .min(1e6);
            (dens * dist, p)
        })
        .collect();
    scored.sort_by(|a, b| b.0.total_cmp(&a.0));
    scored.into_iter().take(k).map(|(_, p)| p).collect()
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
