//! Linear solve for a block's interface Jacobian.
//!
//! A block has one unknown per scalar port side in a loop, and most of them
//! are copies: a spring's end rate is minus a mass's rate, a mode's force is
//! a ratio times the contact's. [`Condensed`] substitutes those away in an
//! order found at sealing and leaves the unknowns that close the loops (the
//! force sum on each mass, the rate sum at a contact) to plain Gaussian
//! elimination with partial pivoting. It is the same linear system, so the
//! Newton iteration is the one a dense solve of every unknown would run.
//!
//! The factors are kept apart from the solve, so a block whose Jacobian does
//! not change from step to step is factored once.

use core::ops::Range;

use crate::wiring::blocks_in_dependency_order;

/// A matrix `I + A` with the pattern of `A` fixed, solved by substitution
/// where the pattern allows it and densely for the rest.
///
/// The unknowns are split at planning. Those that *follow* take places
/// `0..follow` in an order in which each depends only on kept unknowns and on
/// followers before it, so their part of the matrix is unit lower triangular
/// and costs a substitution. The *kept* unknowns take the places after; what
/// remains for them once the followers are substituted is a dense matrix of
/// their own size, factored with partial pivoting.
#[derive(Default)]
pub(crate) struct Condensed {
    /// How many unknowns follow.
    follow: usize,
    /// Place of each unknown.
    place: Vec<usize>,
    /// Per place, where its column of `A` is in `rows` and `values`.
    columns: Vec<Range<usize>>,
    /// Row of each entry, as a place, ascending within a column.
    rows: Vec<usize>,
    values: Vec<f64>,
    /// Per kept unknown, where the followers its column reaches are in
    /// `reached`: their places, ascending.
    reach: Vec<Range<usize>>,
    reached: Vec<usize>,
    /// The kept unknowns' matrix with the followers substituted, row-major;
    /// after [`Self::factor`], its factors.
    dense: Vec<f64>,
    pivots: Vec<usize>,
    /// One value per place; zero between calls.
    work: Vec<f64>,
    /// One value per follower; zero between calls.
    carry: Vec<f64>,
}

impl Condensed {
    /// Plan for a matrix whose column `j` may have entries, besides the
    /// identity's, in the rows `pattern[j]`.
    ///
    /// An unknown can follow if nothing it depends on depends on it. Where
    /// unknowns depend on each other in a ring, the ones that depend on the
    /// most members of the ring are kept, no two of them neighbours in one
    /// round, until no ring is left; of a ring that is one plain cycle, one
    /// is kept. For masses joined by springs that keeps the force sum on
    /// each mass and lets every spring end follow; for modes that share one
    /// contact it keeps the contact's rate alone. The choice is a rule, not
    /// a search for the fewest: keeping more than needed costs time only.
    pub(crate) fn plan(pattern: &[Vec<usize>]) -> Self {
        let m = pattern.len();
        let (follows, order) = followers(pattern);
        let follow = order.len();
        let mut place = vec![0_usize; m];
        for (at, &unknown) in order.iter().enumerate() {
            place[unknown] = at;
        }
        let mut unknown_at = order;
        for unknown in 0..m {
            if !follows[unknown] {
                place[unknown] = unknown_at.len();
                unknown_at.push(unknown);
            }
        }
        let mut columns = Vec::with_capacity(m);
        let mut rows = Vec::new();
        for (at, &unknown) in unknown_at.iter().enumerate() {
            let start = rows.len();
            rows.extend(pattern[unknown].iter().map(|&row| place[row]));
            rows[start..].sort_unstable();
            debug_assert!(at >= follow || rows[start..].iter().all(|&row| row > at));
            columns.push(start..rows.len());
        }
        // The followers each kept column reaches: those it has entries in,
        // and those they have entries in.
        let mut reach = Vec::with_capacity(m - follow);
        let mut reached = Vec::new();
        let mut seen = vec![false; follow];
        for kept in follow..m {
            let start = reached.len();
            let mut next = start;
            let mut visit = |row: usize, reached: &mut Vec<usize>| {
                if row < follow && !seen[row] {
                    seen[row] = true;
                    reached.push(row);
                }
            };
            for &row in &rows[columns[kept].clone()] {
                visit(row, &mut reached);
            }
            while next < reached.len() {
                let at = reached[next];
                next += 1;
                for &row in &rows[columns[at].clone()] {
                    visit(row, &mut reached);
                }
            }
            reached[start..].sort_unstable();
            for &at in &reached[start..] {
                seen[at] = false;
            }
            reach.push(start..reached.len());
        }
        let kept = m - follow;
        Self {
            follow,
            place,
            values: vec![0.0; rows.len()],
            columns,
            rows,
            reach,
            reached,
            dense: vec![0.0; kept * kept],
            pivots: vec![0; kept],
            work: vec![0.0; m],
            carry: vec![0.0; follow],
        }
    }
}

/// Which unknowns of `pattern` follow, and the order to substitute them in:
/// each after everything it depends on. See [`Condensed::plan`].
fn followers(pattern: &[Vec<usize>]) -> (Vec<bool>, Vec<usize>) {
    let m = pattern.len();
    // For each row, the columns with an entry in it: what it depends on.
    let mut feeders: Vec<Vec<usize>> = vec![Vec::new(); m];
    for (column, rows) in pattern.iter().enumerate() {
        for &row in rows {
            feeders[row].push(column);
        }
    }
    let mut follows = vec![true; m];
    loop {
        let dependents: Vec<Vec<usize>> = (0..m)
            .map(|column| {
                if follows[column] {
                    pattern[column].iter().copied().filter(|&row| follows[row]).collect()
                } else {
                    Vec::new()
                }
            })
            .collect();
        let groups = blocks_in_dependency_order(&dependents);
        let mut ring = vec![usize::MAX; m];
        let mut deferred = vec![false; m];
        let mut any = false;
        for (id, group) in groups.iter().enumerate() {
            if group.len() == 1 && !dependents[group[0]].contains(&group[0]) {
                continue;
            }
            any = true;
            for &unknown in group {
                ring[unknown] = id;
            }
            let within =
                |unknown: usize| feeders[unknown].iter().filter(|&&from| ring[from] == id).count();
            let most = group.iter().map(|&unknown| within(unknown)).max().unwrap_or(0);
            for &unknown in group {
                if deferred[unknown] || within(unknown) != most {
                    continue;
                }
                follows[unknown] = false;
                if most == 1 {
                    // Every member depends on one other: a plain cycle,
                    // open once any one of it is kept.
                    break;
                }
                for &near in feeders[unknown].iter().chain(&pattern[unknown]) {
                    if ring[near] == id {
                        deferred[near] = true;
                    }
                }
            }
        }
        if !any {
            let order = groups
                .into_iter()
                .map(|group| group[0])
                .filter(|&unknown| follows[unknown])
                .collect();
            return (follows, order);
        }
    }
}

impl Condensed {
    /// How many unknowns are solved together; the rest follow from them.
    pub(crate) fn kept(&self) -> usize {
        self.place.len() - self.follow
    }

    /// Where the entry of `row` in `column` is kept, for [`Self::subtract`].
    /// The entry must be in the planned pattern.
    pub(crate) fn slot(&self, row: usize, column: usize) -> usize {
        let range = self.columns[self.place[column]].clone();
        let at = self.rows[range.clone()]
            .binary_search(&self.place[row])
            .expect("an entry of the planned pattern");
        range.start + at
    }

    /// Make the matrix the identity again.
    pub(crate) fn clear(&mut self) {
        self.values.fill(0.0);
    }

    /// Subtract `value` from the entry at `slot`.
    pub(crate) fn subtract(&mut self, slot: usize, value: f64) {
        self.values[slot] -= value;
    }

    /// Substitute the followers out of the kept unknowns' equations and
    /// factor what remains. Returns `false` if that matrix is singular; the
    /// factors must then not be used.
    pub(crate) fn factor(&mut self) -> bool {
        let follow = self.follow;
        let kept = self.kept();
        self.dense.fill(0.0);
        for index in 0..kept {
            self.dense[index * kept + index] = 1.0;
            // This column, then what the followers it reaches pass on, each
            // once every follower before it has.
            for entry in self.columns[follow + index].clone() {
                self.work[self.rows[entry]] += self.values[entry];
            }
            for &at in &self.reached[self.reach[index].clone()] {
                let through = core::mem::take(&mut self.work[at]);
                if through != 0.0 {
                    for entry in self.columns[at].clone() {
                        self.work[self.rows[entry]] -= self.values[entry] * through;
                    }
                }
            }
            for row in 0..kept {
                self.dense[row * kept + index] += core::mem::take(&mut self.work[follow + row]);
            }
        }
        factor_dense(&mut self.dense, &mut self.pivots, kept)
    }

    /// Solve `(I + A) x = b` with the factors of [`Self::factor`]; `b`
    /// becomes `x`. Returns `false` if the solution is not finite.
    pub(crate) fn solve(&mut self, b: &mut [f64]) -> bool {
        let follow = self.follow;
        let kept = self.kept();
        debug_assert_eq!(b.len(), self.place.len());
        for (&value, &at) in b.iter().zip(&self.place) {
            self.work[at] = value;
        }
        // Followers in order, as if the kept unknowns were zero; each passes
        // its value on to the equations after it.
        for at in 0..follow {
            let value = self.work[at];
            if value != 0.0 {
                for entry in self.columns[at].clone() {
                    self.work[self.rows[entry]] -= self.values[entry] * value;
                }
            }
        }
        let solved = solve_factored(&self.dense, &self.pivots, &mut self.work[follow..], kept);
        // What the kept unknowns add to each follower, by the same
        // substitution.
        for index in 0..kept {
            let value = self.work[follow + index];
            if value != 0.0 {
                for entry in self.columns[follow + index].clone() {
                    let row = self.rows[entry];
                    if row >= follow {
                        break;
                    }
                    self.carry[row] += self.values[entry] * value;
                }
            }
        }
        for at in 0..follow {
            let value = core::mem::take(&mut self.carry[at]);
            self.work[at] -= value;
            if value != 0.0 {
                for entry in self.columns[at].clone() {
                    let row = self.rows[entry];
                    if row >= follow {
                        break;
                    }
                    self.carry[row] -= self.values[entry] * value;
                }
            }
        }
        let mut finite = solved;
        for (value, &at) in b.iter_mut().zip(&self.place) {
            *value = core::mem::take(&mut self.work[at]);
            finite &= value.is_finite();
        }
        finite
    }
}

/// Factor a row-major `n x n` matrix in place as `P A = L U` with partial
/// pivoting. `a` then holds `U` on and above the diagonal and the multipliers
/// of the unit lower triangle `L` below it; `pivots[k]` is the row exchanged
/// with row `k`.
///
/// Returns `false` if a pivot is zero or not finite; `a` is then partly
/// eliminated and must not be used.
pub(crate) fn factor_dense(a: &mut [f64], pivots: &mut [usize], n: usize) -> bool {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(pivots.len(), n);
    for column in 0..n {
        let mut pivot_row = column;
        let mut pivot = a[column * n + column].abs();
        for row in column + 1..n {
            let candidate = a[row * n + column].abs();
            if candidate > pivot {
                pivot = candidate;
                pivot_row = row;
            }
        }
        if !(pivot.is_finite() && pivot > 0.0) {
            return false;
        }
        pivots[column] = pivot_row;
        if pivot_row != column {
            // Whole rows, so the multipliers already stored move with them.
            for k in 0..n {
                a.swap(column * n + k, pivot_row * n + k);
            }
        }
        let diagonal = a[column * n + column];
        for row in column + 1..n {
            let factor = a[row * n + column] / diagonal;
            a[row * n + column] = factor;
            if factor != 0.0 {
                for k in column + 1..n {
                    a[row * n + k] -= factor * a[column * n + k];
                }
            }
        }
    }
    true
}

/// Solve `A x = b` with the factors of [`factor_dense`]; `b` becomes `x`.
///
/// Returns `false` if the solution is not finite.
pub(crate) fn solve_factored(a: &[f64], pivots: &[usize], b: &mut [f64], n: usize) -> bool {
    debug_assert_eq!(a.len(), n * n);
    debug_assert_eq!(pivots.len(), n);
    debug_assert_eq!(b.len(), n);
    for (row, &exchanged) in pivots.iter().enumerate() {
        b.swap(row, exchanged);
    }
    for column in 0..n {
        for row in column + 1..n {
            let factor = a[row * n + column];
            if factor != 0.0 {
                b[row] -= factor * b[column];
            }
        }
    }
    for row in (0..n).rev() {
        let mut sum = b[row];
        for k in row + 1..n {
            sum -= a[row * n + k] * b[k];
        }
        b[row] = sum / a[row * n + row];
    }
    b.iter().all(|value| value.is_finite())
}

#[cfg(test)]
mod tests {
    use super::{Condensed, factor_dense, solve_factored};

    /// Deterministic values in `(-0.5, 0.5)`.
    fn values(mut seed: u64) -> impl FnMut() -> f64 {
        move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        }
    }

    /// `I + A` on `pattern` with entries of size up to `size`, planned and
    /// as a dense row-major matrix.
    fn filled(pattern: &[Vec<usize>], size: f64, next: &mut impl FnMut() -> f64) -> (Condensed, Vec<f64>) {
        let m = pattern.len();
        let mut condensed = Condensed::plan(pattern);
        let mut dense = vec![0.0; m * m];
        for diagonal in 0..m {
            dense[diagonal * m + diagonal] = 1.0;
        }
        for (column, rows) in pattern.iter().enumerate() {
            for &row in rows {
                let value = 2.0 * size * next();
                condensed.subtract(condensed.slot(row, column), value);
                dense[row * m + column] -= value;
            }
        }
        (condensed, dense)
    }

    /// Masses in a line with a spring between neighbours, as the kernel
    /// sees them: unknown `i < n` is the force sum on mass `i`; spring `s`
    /// has an end rate at each of its two masses. A spring's forces depend
    /// on both its ends, an end on its mass's force.
    fn line(n: usize) -> Vec<Vec<usize>> {
        let mut pattern = vec![Vec::new(); n + 2 * (n - 1)];
        for spring in 0..n - 1 {
            for (end, mass) in [(n + 2 * spring, spring), (n + 2 * spring + 1, spring + 1)] {
                pattern[mass].push(end);
                pattern[end] = vec![spring, spring + 1];
            }
        }
        pattern
    }

    /// One hub that sums `n` spokes, each of which copies the hub.
    fn star(n: usize) -> Vec<Vec<usize>> {
        let mut pattern = vec![(1..=n).collect::<Vec<_>>()];
        pattern.extend((0..n).map(|_| vec![0]));
        pattern
    }

    #[test]
    fn the_unknowns_kept_are_the_ones_that_close_the_loops() {
        // Each mass's force sum, and none of the 2 (n - 1) spring ends.
        for n in [2_usize, 3, 8, 48] {
            let plan = Condensed::plan(&line(n));
            assert_eq!(plan.kept(), n, "line of {n}");
            assert!((0..n).all(|mass| plan.place[mass] >= plan.follow));
        }
        // The hub's sum alone.
        let plan = Condensed::plan(&star(200));
        assert_eq!(plan.kept(), 1);
        assert_eq!(plan.place[0], 200);
        // A chain has no loop, so nothing is kept; a plain cycle is opened
        // at one unknown, and an unknown that depends on itself is kept.
        let chain: Vec<Vec<usize>> = (0..6).map(|at| if at < 5 { vec![at + 1] } else { vec![] }).collect();
        assert_eq!(Condensed::plan(&chain).kept(), 0);
        let cycle: Vec<Vec<usize>> = (0..7).map(|at| vec![(at + 1) % 7]).collect();
        assert_eq!(Condensed::plan(&cycle).kept(), 1);
        let own = vec![vec![0, 1], vec![]];
        let plan = Condensed::plan(&own);
        assert_eq!((plan.kept(), plan.place[0]), (1, 1));
        // The same pattern gives the same plan.
        let (a, b) = (Condensed::plan(&line(9)), Condensed::plan(&line(9)));
        assert_eq!((a.place, a.rows, a.reached), (b.place, b.rows, b.reached));
        assert_eq!(Condensed::plan(&[]).kept(), 0);
    }

    #[test]
    fn the_rule_keeps_no_two_neighbours_in_a_round_and_counts_within_a_ring() {
        // Every unknown depends on every other. One is kept and the rest,
        // its neighbours, wait for the next round; in the end one follows.
        // Keeping every unknown of the largest count at once would keep all.
        for m in [3_usize, 6] {
            let full: Vec<Vec<usize>> =
                (0..m).map(|column| (0..m).filter(|&row| row != column).collect()).collect();
            assert_eq!(Condensed::plan(&full).kept(), m - 1, "{m} unknowns");
        }
        // A pair that depend on each other, and a cycle of four of which
        // two also depend on the pair. The cycle is a plain cycle whatever
        // feeds it from outside, so one of it is kept: two in all. Counting
        // what an unknown depends on outside its own ring would keep both
        // of the fed ones.
        let fed = vec![vec![1, 2, 4], vec![0, 2, 4], vec![3], vec![4], vec![5], vec![2]];
        let plan = Condensed::plan(&fed);
        assert_eq!(plan.kept(), 2);
        assert!(plan.place[0] >= plan.follow && plan.place[2] >= plan.follow);
    }

    #[test]
    fn patterns_drawn_at_random_solve_what_dense_elimination_solves() {
        let mut next = values(0x1234_5678_9ABC_DEF1);
        let mut kept_some = 0_u32;
        let mut followed_some = 0_u32;
        for trial in 0..600_usize {
            let m = 1 + trial % 14;
            let pattern: Vec<Vec<usize>> = if trial % 2 == 0 {
                // Any entry, about one in four.
                (0..m).map(|_| (0..m).filter(|_| next() > 0.25).collect()).collect()
            } else {
                // As the kernel makes them: each unknown is an input of a
                // component, each row reads one or two components, and a
                // column has an entry in every row that reads its owner.
                let components = 1 + trial % 5;
                let pick = |value: f64| ((value + 0.5) * components as f64) as usize % components;
                let owner: Vec<usize> = (0..m).map(|_| pick(next())).collect();
                let reads: Vec<[usize; 2]> = (0..m).map(|_| [pick(next()), pick(next())]).collect();
                (0..m)
                    .map(|column| (0..m).filter(|&row| reads[row].contains(&owner[column])).collect())
                    .collect()
            };
            let (mut condensed, matrix) = filled(&pattern, 0.1, &mut next);
            kept_some += u32::from(condensed.kept() > 0);
            followed_some += u32::from(condensed.follow > 0);
            let mut factors = matrix.clone();
            let mut pivots = vec![0; m];
            assert!(factor_dense(&mut factors, &mut pivots, m), "trial {trial}");
            assert!(condensed.factor(), "trial {trial}");
            let right: Vec<f64> = (0..m).map(|_| next()).collect();
            let mut got = right.clone();
            assert!(condensed.solve(&mut got));
            let mut want = right;
            assert!(solve_factored(&factors, &pivots, &mut want, m));
            for (a, b) in got.iter().zip(&want) {
                assert!((a - b).abs() <= 1.0e-12 * (1.0 + b.abs()), "trial {trial}: {a} vs {b}");
            }
        }
        // Both parts of the solve were exercised.
        assert!(kept_some > 200 && followed_some > 200, "{kept_some} {followed_some}");
    }

    #[test]
    fn substitution_and_the_kept_factors_solve_what_dense_elimination_solves() {
        let mut next = values(0x2545_F491_4F6C_DD1D);
        // Every unknown depends on every other: nearly all are kept.
        let full: Vec<Vec<usize>> =
            (0..6).map(|column| (0..6).filter(|&row| row != column).collect()).collect();
        // Rings that share members, a tail that only follows, self entries.
        let mixed = vec![
            vec![1, 4], vec![2, 0], vec![0, 2, 5], vec![3, 0], vec![3], vec![6], vec![7], vec![],
        ];
        let cycle: Vec<Vec<usize>> = (0..7).map(|at| vec![(at + 1) % 7]).collect();
        let chain: Vec<Vec<usize>> = (0..6).map(|at| if at < 5 { vec![at + 1] } else { vec![] }).collect();
        // One kept unknown that reaches a follower both directly and
        // through two others: what they pass on has to arrive before that
        // follower passes its own value on.
        let rejoined = vec![vec![1, 3, 4], vec![2], vec![3], vec![0], vec![0]];
        assert_eq!(Condensed::plan(&rejoined).kept(), 1);
        for (pattern, size) in [
            (line(2), 0.5), (line(12), 0.5), (star(30), 0.15), (full, 0.15), (mixed, 0.4),
            (cycle, 0.5), (chain, 0.5), (rejoined, 0.5), (Vec::new(), 0.5),
        ] {
            let m = pattern.len();
            let (mut condensed, matrix) = filled(&pattern, size, &mut next);
            assert!(condensed.factor());
            let mut factors = matrix.clone();
            let mut pivots = vec![0; m];
            assert!(factor_dense(&mut factors, &mut pivots, m));
            for _ in 0..3 {
                let right: Vec<f64> = (0..m).map(|_| next()).collect();
                let mut got = right.clone();
                assert!(condensed.solve(&mut got));
                let mut want = right.clone();
                assert!(solve_factored(&factors, &pivots, &mut want, m));
                for (a, b) in got.iter().zip(&want) {
                    assert!((a - b).abs() <= 1.0e-13 * (1.0 + b.abs()), "{a} vs {b} ({m} unknowns)");
                }
                // The factors are not used up: the same right-hand side
                // again gives the same bits, and the workspaces are clean.
                let mut again = right.clone();
                assert!(condensed.solve(&mut again));
                assert_eq!(
                    again.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    got.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
                assert!(condensed.work.iter().chain(&condensed.carry).all(|&x| x == 0.0));
            }
            // Cleared and filled with other values, it is another matrix.
            condensed.clear();
            assert!(condensed.factor());
            let mut same: Vec<f64> = (0..m).map(|_| next()).collect();
            let before = same.clone();
            assert!(condensed.solve(&mut same));
            assert_eq!(same, before);
        }
    }

    #[test]
    fn a_singular_kept_matrix_and_an_overflowing_solution_are_refused() {
        // x0 = x1 and x1 = x0: one is kept, and its equation is 0 = b.
        let pattern = vec![vec![1], vec![0]];
        let mut condensed = Condensed::plan(&pattern);
        condensed.subtract(condensed.slot(1, 0), 1.0);
        condensed.subtract(condensed.slot(0, 1), 1.0);
        assert!(!condensed.factor());
        // A follower that overflows.
        let chain = vec![vec![1], vec![]];
        let mut condensed = Condensed::plan(&chain);
        condensed.subtract(condensed.slot(1, 0), f64::MAX);
        assert!(condensed.factor());
        let mut b = vec![4.0, f64::MAX];
        assert!(!condensed.solve(&mut b));
        // It can be used again afterwards.
        let mut b = vec![0.0, 3.0];
        assert!(condensed.solve(&mut b));
        assert_eq!(b, [0.0, 3.0]);

        // A plan that refused a singular matrix, or was factored and never
        // solved, is as good as new once cleared and filled again.
        let pattern = line(5);
        let mut next = values(0x0DDB_A11C_0FFE_E123);
        let (mut fresh, _) = filled(&pattern, 0.4, &mut next);
        let mut next = values(0x0DDB_A11C_0FFE_E123);
        let mut used = Condensed::plan(&pattern);
        for (column, rows) in pattern.iter().enumerate() {
            for &row in rows {
                // Entries that are not numbers reach every workspace.
                used.subtract(used.slot(row, column), f64::NAN);
            }
        }
        assert!(!used.factor());
        for _ in 0..2 {
            used.clear();
            for (column, rows) in pattern.iter().enumerate() {
                for &row in rows {
                    used.subtract(used.slot(row, column), 0.8 * next());
                }
            }
            assert!(used.factor());
            next = values(0x0DDB_A11C_0FFE_E123);
        }
        assert!(fresh.factor());
        let right: Vec<f64> = (0..pattern.len()).map(|at| 0.3 + at as f64).collect();
        let (mut a, mut b) = (right.clone(), right);
        assert!(used.solve(&mut a) && fresh.solve(&mut b));
        assert_eq!(
            a.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
        );
    }

    /// Factor and solve in one call; `a` is replaced by its factors.
    fn solve_dense(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
        let mut pivots = vec![0; n];
        factor_dense(a, &mut pivots, n) && solve_factored(a, &pivots, b, n)
    }

    /// Elimination applied to the right-hand side as it goes, with no kept
    /// factors: the solve this module had before the factors were kept.
    fn eliminate(a: &mut [f64], b: &mut [f64], n: usize) -> bool {
        for column in 0..n {
            let mut pivot_row = column;
            let mut pivot = a[column * n + column].abs();
            for row in column + 1..n {
                let candidate = a[row * n + column].abs();
                if candidate > pivot {
                    pivot = candidate;
                    pivot_row = row;
                }
            }
            if !(pivot.is_finite() && pivot > 0.0) {
                return false;
            }
            if pivot_row != column {
                for k in 0..n {
                    a.swap(column * n + k, pivot_row * n + k);
                }
                b.swap(column, pivot_row);
            }
            let diagonal = a[column * n + column];
            for row in column + 1..n {
                let factor = a[row * n + column] / diagonal;
                if factor != 0.0 {
                    for k in column..n {
                        a[row * n + k] -= factor * a[column * n + k];
                    }
                    b[row] -= factor * b[column];
                }
            }
        }
        for row in (0..n).rev() {
            let mut sum = b[row];
            for k in row + 1..n {
                sum -= a[row * n + k] * b[k];
            }
            b[row] = sum / a[row * n + row];
        }
        b.iter().all(|value| value.is_finite())
    }

    #[test]
    fn solves_a_system_that_needs_pivoting() {
        // First pivot is zero without a row exchange.
        let mut a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 0.0, 3.0, 0.0, 1.0];
        // b = a * [1, 2, 1].
        let mut b = vec![5.0, 3.0, 4.0];
        assert!(solve_dense(&mut a, &mut b, 3));
        let expected = [1.0, 2.0, 1.0];
        for (got, want) in b.iter().zip(expected) {
            assert!((got - want).abs() <= 1.0e-14, "{got} vs {want}");
        }
    }

    #[test]
    fn refuses_a_singular_matrix() {
        let mut a = vec![1.0, 2.0, 2.0, 4.0];
        let mut b = vec![1.0, 2.0];
        assert!(!solve_dense(&mut a, &mut b, 2));
    }

    #[test]
    fn an_empty_system_is_solved() {
        assert!(solve_dense(&mut [], &mut [], 0));
    }

    #[test]
    fn kept_factors_solve_every_right_hand_side_to_the_bit_of_one_elimination() {
        // Matrices that exchange rows at several columns, with zeros that
        // make some multipliers vanish. Deterministic pseudo-random entries.
        let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
        };
        for n in [1_usize, 2, 5, 12] {
            let matrix: Vec<f64> = (0..n * n)
                .map(|entry| {
                    let value = next();
                    // A seventh of the entries are zero (never a whole
                    // column at these sizes), the diagonal small.
                    if entry % 7 == 3 {
                        0.0
                    } else if entry % (n + 1) == 0 {
                        1.0e-3 * value
                    } else {
                        value
                    }
                })
                .collect();
            let mut factors = matrix.clone();
            let mut pivots = vec![0; n];
            assert!(factor_dense(&mut factors, &mut pivots, n));
            assert!(n < 5 || pivots.iter().enumerate().any(|(row, &to)| row != to));
            for _ in 0..3 {
                let right: Vec<f64> = (0..n).map(|_| next()).collect();
                let mut kept = right.clone();
                assert!(solve_factored(&factors, &pivots, &mut kept, n));
                let mut once = right.clone();
                assert!(eliminate(&mut matrix.clone(), &mut once, n));
                for (a, b) in kept.iter().zip(&once) {
                    assert_eq!(a.to_bits(), b.to_bits());
                }
                // And it is a solution.
                for row in 0..n {
                    let product: f64 = (0..n).map(|k| matrix[row * n + k] * kept[k]).sum();
                    assert!((product - right[row]).abs() <= 1.0e-9, "{product} vs {}", right[row]);
                }
            }
        }
    }
}
