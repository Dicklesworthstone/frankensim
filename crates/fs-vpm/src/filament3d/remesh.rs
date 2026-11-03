//! Conservative geometric refinement of the connected three-dimensional sheet.

use super::{FilamentWake, Refusal, ShedRow};

/// Result of one bounded streamwise refinement; no physical time has elapsed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StreamwiseRefinement {
    /// Retained row count before refinement.
    pub rows_before: usize,
    /// Retained row count after refinement.
    pub rows_after: usize,
    /// Longest streamwise edge before refinement [m], or zero for no rows.
    pub longest_edge_before_m: f64,
    /// Longest streamwise edge after refinement [m], or zero for no rows.
    pub longest_edge_after_m: f64,
}

impl FilamentWake {
    /// Subdivide every streamwise row interval whose longest edge exceeds
    /// `max_edge_length_m`, inserting complete spanwise node lines. Each
    /// interval uses uniform subdivisions, including its zero-strength edges.
    /// The returned longest edge meets the target within 64 machine epsilons.
    ///
    /// Each child quad retains its parent's circulation, base core and age.
    /// New interior spanwise strengths therefore cancel exactly; existing
    /// straight trailing filaments are only subdivided. The physical velocity,
    /// linear impulse and angular impulse remain unchanged to roundoff. Ages
    /// must not be interpolated: doing so would change the finite-core field.
    /// Existing node lines, the lifting line, ticks and an attached far-field
    /// boundary remain fixed. This operation does not advance time, coarsen,
    /// or model viscous diffusion or previously unresolved curvature.
    ///
    /// Work is linear in original plus inserted lattice size. Both the caller's
    /// `max_new_rows` and the certificate row cap must admit the entire change.
    /// Cancellation is polled before admission, for each original/inserted row,
    /// and before publication. Every refusal leaves the wake unchanged.
    ///
    /// # Errors
    /// Invalid geometry, target, age/core or interface; insufficient row budget;
    /// allocation failure; subdivision below coordinate resolution; cancellation.
    pub fn refine_streamwise(
        &mut self,
        max_edge_length_m: f64,
        max_new_rows: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StreamwiseRefinement, Refusal> {
        checkpoint(&mut cancelled)?;
        self.cert.admit()?;
        if !max_edge_length_m.is_finite() || max_edge_length_m <= 0.0 {
            return Err(refusal(
                "wake-refinement-invalid",
                "edge target must be finite and positive",
            ));
        }
        let available = self
            .cert
            .max_rows
            .checked_sub(self.rows.len())
            .ok_or_else(|| {
                refusal(
                    "wake-refinement-budget",
                    "existing rows exceed the certificate",
                )
            })?
            .min(max_new_rows);
        self.admit_geometry()?;
        if self.downstream_gamma.as_ref().is_some_and(|gamma| {
            gamma.len() != self.cert.n_stations || gamma.iter().any(|g| !g.is_finite())
        }) {
            return Err(refusal(
                "wake-refinement-invalid",
                "invalid attached interface circulation",
            ));
        }

        let mut divisions = Vec::new();
        reserve(&mut divisions, self.rows.len())?;
        let mut added = 0;
        let mut longest_before: f64 = 0.0;
        let mut longest_after: f64 = 0.0;
        for (r, row) in self.rows.iter().enumerate() {
            checkpoint(&mut cancelled)?;
            if !row.core2_m2.is_finite()
                || row.core2_m2 < 0.0
                || !row.age_ticks.is_finite()
                || row.age_ticks < 0.0
            {
                return Err(refusal(
                    "wake-refinement-invalid",
                    "retained age and core must be finite and nonnegative",
                ));
            }
            let next = self
                .rows
                .get(r + 1)
                .map_or(&self.line_nodes, |row| &row.nodes);
            let length = longest_edge(&row.nodes, next)?;
            longest_before = longest_before.max(length);
            let ratio = length / max_edge_length_m;
            let pieces = if within_target(ratio) {
                1.0
            } else {
                ratio.ceil()
            };
            // Compare while still floating point: an extreme ratio must not
            // saturate an integer cast or overflow an allocation count.
            if pieces > (available - added + 1) as f64 {
                return Err(refusal(
                    "wake-refinement-budget",
                    "complete refinement exceeds the caller or certificate row cap",
                ));
            }
            let pieces = pieces as usize;
            added += pieces - 1;
            divisions.push(pieces);
            if pieces == 1 {
                longest_after = longest_after.max(length);
            }
        }
        let mut report = StreamwiseRefinement {
            rows_before: self.rows.len(),
            rows_after: self.rows.len() + added,
            longest_edge_before_m: longest_before,
            longest_edge_after_m: longest_after,
        };
        if added == 0 {
            checkpoint(&mut cancelled)?;
            return Ok(report);
        }

        // Prepare only inserted buffers. The original rows are moved intact
        // after the last fallible operation, without cloning the existing mesh.
        let mut inserted: Vec<ShedRow> = Vec::new();
        reserve(&mut inserted, added)?;
        let mut output = Vec::new();
        reserve(&mut output, report.rows_after)?;
        for (r, (&pieces, row)) in divisions.iter().zip(&self.rows).enumerate() {
            if pieces == 1 {
                continue;
            }
            let next = self
                .rows
                .get(r + 1)
                .map_or(&self.line_nodes, |row| &row.nodes);
            for part in 1..pieces {
                checkpoint(&mut cancelled)?;
                let fraction = part as f64 / pieces as f64;
                let mut nodes = Vec::new();
                reserve(&mut nodes, row.nodes.len())?;
                for (a, b) in row.nodes.iter().zip(next) {
                    nodes.push(std::array::from_fn(|axis| {
                        a[axis] + fraction * (b[axis] - a[axis])
                    }));
                }
                let previous = if part == 1 {
                    &row.nodes
                } else {
                    &inserted.last().expect("previous subdivision exists").nodes
                };
                let length = longest_edge(previous, &nodes)?;
                admit_length(length, max_edge_length_m)?;
                longest_after = longest_after.max(length);
                let mut gamma = Vec::new();
                reserve(&mut gamma, row.gamma.len())?;
                gamma.extend_from_slice(&row.gamma);
                inserted.push(ShedRow {
                    nodes,
                    gamma,
                    core2_m2: row.core2_m2,
                    age_ticks: row.age_ticks,
                });
            }
            let length = longest_edge(&inserted.last().expect("split inserted a row").nodes, next)?;
            admit_length(length, max_edge_length_m)?;
            longest_after = longest_after.max(length);
        }
        report.longest_edge_after_m = longest_after;
        checkpoint(&mut cancelled)?;
        let mut inserted = inserted.into_iter();
        for (row, pieces) in std::mem::take(&mut self.rows).into_iter().zip(divisions) {
            output.push(row);
            for _ in 1..pieces {
                output.push(inserted.next().expect("all subdivisions prepared"));
            }
        }
        self.rows = output;
        Ok(report)
    }
}

fn longest_edge(a: &[[f64; 3]], b: &[[f64; 3]]) -> Result<f64, Refusal> {
    let mut longest: f64 = 0.0;
    for (a, b) in a.iter().zip(b) {
        let length = (b[0] - a[0]).hypot(b[1] - a[1]).hypot(b[2] - a[2]);
        if !length.is_finite() {
            return Err(refusal(
                "wake-refinement-unrepresentable",
                "streamwise edge exceeds the finite coordinate domain",
            ));
        }
        longest = longest.max(length);
    }
    Ok(longest)
}

fn admit_length(length: f64, target: f64) -> Result<(), Refusal> {
    if !within_target(length / target) {
        return Err(refusal(
            "wake-refinement-unrepresentable",
            "requested edge length is below coordinate resolution",
        ));
    }
    Ok(())
}

fn within_target(ratio: f64) -> bool {
    ratio <= 1.0 + 64.0 * f64::EPSILON
}

fn reserve<T>(items: &mut Vec<T>, count: usize) -> Result<(), Refusal> {
    items.try_reserve_exact(count).map_err(|_| {
        refusal(
            "wake-refinement-allocation",
            "refinement scratch allocation failed",
        )
    })
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "wake-refinement-cancelled",
            "streamwise refinement was cancelled",
        ))
    } else {
        Ok(())
    }
}

fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use a finite wake, a resolvable positive edge target and sufficient admitted rows"
                .into(),
        ],
    }
}
