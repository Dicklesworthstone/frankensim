//! Bracket one monotone feedback coordinate after substituting its graph.
//! The ordinary convergence and conservation checks still decide acceptance.

use super::*;

pub(super) struct ScalarPlan {
    /// Position in the block's unknowns of the one input left open.
    row: usize,
    /// Removing that input's incoming dependencies makes this order acyclic.
    order: Vec<usize>,
}

impl SystemBuilder {
    pub(super) fn scalar_plan(
        &self,
        block: &Block,
        drivers: &[Driver],
        input_owner: &[usize],
        output_owner: &[usize],
    ) -> Option<ScalarPlan> {
        if block
            .nodes
            .iter()
            .any(|&node| self.nodes[node].class() == Class::NonSmooth)
        {
            return None;
        }
        for &node in &block.nodes {
            if self.nodes[node].class() != Class::Monotone1d || self.nodes[node].has_tangent() {
                continue;
            }
            let mut own = block
                .unknowns
                .iter()
                .enumerate()
                .filter(|&(_, &input)| input_owner[input] == node);
            let Some((row, &tear)) = own.next() else {
                continue;
            };
            if own.next().is_some() {
                continue;
            }
            let mut edges = vec![Vec::new(); block.nodes.len()];
            let mut sources = Vec::new();
            for &input in &block.unknowns {
                if input == tear {
                    continue;
                }
                let to = block.nodes.binary_search(&input_owner[input]).ok()?;
                sources.clear();
                drivers[input].sources(&mut sources);
                for &output in &sources {
                    if let Ok(from) = block.nodes.binary_search(&output_owner[output]) {
                        if !edges[from].contains(&to) {
                            edges[from].push(to);
                        }
                    }
                }
            }
            let groups = blocks_in_dependency_order(&edges);
            if groups
                .iter()
                .any(|group| group.len() != 1 || edges[group[0]].contains(&group[0]))
            {
                continue;
            }
            return Some(ScalarPlan {
                row,
                order: groups.iter().map(|group| block.nodes[group[0]]).collect(),
            });
        }
        None
    }
}

#[derive(Clone, Copy)]
struct Point {
    x: f64,
    residual: f64,
}

struct SearchEnd {
    valid: Point,
    // A refused input is a search boundary, never a signed residual.
    limit: Option<f64>,
    distance: f64,
    direction: f64,
    exhausted: bool,
}

fn midpoint(a: f64, b: f64) -> f64 {
    // Avoid overflow in b-a for large opposite values.
    if a.is_sign_positive() == b.is_sign_positive() {
        a + 0.5 * (b - a)
    } else {
        0.5 * a + 0.5 * b
    }
}

impl System {
    /// Evaluate the actual graph with one prescribed feedback input. Every
    /// other loop input follows its driver exactly, in the sealed order.
    fn scalar_value(
        &mut self,
        plan: &ScalarPlan,
        unknowns: &[usize],
        x: f64,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<Point, f64>, SystemError> {
        cx.checkpoint().map_err(|_| SystemError::Cancelled)?;
        let tear = unknowns[plan.row];
        self.inputs[tear] = x;
        for &node in &plan.order {
            for &input in unknowns {
                if input != tear && self.input_owner[input] == node {
                    self.inputs[input] =
                        self.drivers[input].value(&self.outputs, &self.boundary_values);
                }
            }
            if let Some(shorter) = self.propose(node, dt_s, cx)? {
                return Ok(Err(shorter));
            }
        }
        let target = self.drivers[tear].value(&self.outputs, &self.boundary_values);
        Ok(Ok(Point {
            x,
            residual: x - target,
        }))
    }

    // A seed criterion only. The block solver subsequently checks every wire
    // and carries its residual through locally evaluated slopes.
    fn scalar_near(&self, point: Point) -> bool {
        point.residual.is_finite()
            && point.residual.abs()
                <= self.config.absolute_tolerance
                    + self.config.relative_tolerance
                        * point.x.abs().max((point.x - point.residual).abs())
    }

    /// Recover physical-domain refusals of an exploratory probe. Invalid
    /// contracts and cancellation still stop the complete trial.
    fn scalar_probe(
        &mut self,
        plan: &ScalarPlan,
        unknowns: &[usize],
        x: f64,
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Option<Point>, SystemError> {
        match self.scalar_value(plan, unknowns, x, dt_s, cx) {
            Ok(Ok(point)) if point.residual.is_finite() => Ok(Some(point)),
            Ok(_)
            | Err(SystemError::Component {
                phase: Phase::Propose,
                ..
            }) => Ok(None),
            Err(error) => Err(error),
        }
    }

    pub(super) fn scalar_seed(
        &mut self,
        plan: &ScalarPlan,
        unknowns: &[usize],
        nodes: &[usize],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Result<u32, f64>, SystemError> {
        let initial =
            match self.scalar_value(plan, unknowns, self.inputs[unknowns[plan.row]], dt_s, cx)? {
                Ok(point) => point,
                Err(shorter) => return Ok(Err(shorter)),
            };
        let failure = |system: &Self, point: Point, iterations| SystemError::NotConverged {
            nodes: system.names_of(nodes),
            iterations,
            worst: system.input_label[unknowns[plan.row]].clone(),
            residual: point.residual,
        };
        if !initial.residual.is_finite() {
            return Err(failure(self, initial, 0));
        }
        if self.scalar_near(initial) {
            return Ok(Ok(0));
        }
        // A fixed-point correction has this input's own units. A unit-sized
        // floor can put every probe outside a small source-card support.
        // Only round a direction's radius up to one ULP when the correction
        // would otherwise leave that direction at the initial float.
        let distance = initial.residual.abs();
        let direction = -initial.residual.signum();
        let mut ends = [direction, -direction].map(|direction| {
            let distance = if initial.x + direction * distance == initial.x {
                let adjacent = if direction > 0.0 {
                    initial.x.next_up()
                } else {
                    initial.x.next_down()
                };
                (adjacent - initial.x).abs()
            } else {
                distance
            };
            SearchEnd {
                valid: initial,
                limit: None,
                distance,
                direction,
                exhausted: false,
            }
        });
        let mut iterations = 0;
        let mut bracket = None;
        while iterations < self.config.max_iterations && bracket.is_none() {
            let side = iterations as usize % 2;
            if ends.iter().all(|end| end.exhausted) {
                break;
            }
            let end = &mut ends[side];
            iterations += 1;
            if end.exhausted {
                continue;
            }
            let next = match end.limit {
                Some(limit) => midpoint(end.valid.x, limit),
                None => initial.x + end.direction * end.distance,
            };
            if !next.is_finite() || next == end.valid.x {
                end.exhausted = true;
                continue;
            }
            let Some(point) = self.scalar_probe(plan, unknowns, next, dt_s, cx)? else {
                end.limit = Some(next);
                continue;
            };
            if self.scalar_near(point) {
                return Ok(Ok(iterations));
            }
            if point.residual.is_sign_positive() != end.valid.residual.is_sign_positive() {
                bracket = Some((end.valid, point));
            } else {
                end.valid = point;
                end.distance *= 2.0;
            }
        }
        let Some((mut a, mut b)) = bracket else {
            return Err(failure(self, initial, iterations));
        };
        // Illinois regula falsi halves the weight of a repeatedly retained
        // endpoint. Every third refinement bisects to prevent stagnation.
        let (mut wa, mut wb) = (a.residual, b.residual);
        let mut retained = 0_i8;
        let mut refinement = 0_u32;
        while iterations < self.config.max_iterations {
            refinement += 1;
            let magnitude = wa.abs().max(wb.abs());
            let fraction =
                (wa / magnitude).abs() / ((wa / magnitude).abs() + (wb / magnitude).abs());
            let interpolated = (1.0 - fraction) * a.x + fraction * b.x;
            let next = if refinement % 3 == 0
                || !interpolated.is_finite()
                || interpolated == a.x
                || interpolated == b.x
            {
                midpoint(a.x, b.x)
            } else {
                interpolated
            };
            if next == a.x || next == b.x {
                break;
            }
            iterations += 1;
            // A hole inside admitted endpoints is not bridged with a fake
            // residual sign; the owner's refusal remains authoritative.
            let point = match self.scalar_value(plan, unknowns, next, dt_s, cx)? {
                Ok(point) => point,
                Err(shorter) => return Ok(Err(shorter)),
            };
            if !point.residual.is_finite() {
                return Err(failure(self, point, iterations));
            }
            if self.scalar_near(point) {
                return Ok(Ok(iterations));
            }
            if point.residual.is_sign_positive() == a.residual.is_sign_positive() {
                a = point;
                wa = point.residual;
                wb = if retained == 1 { 0.5 * wb } else { b.residual };
                retained = 1;
            } else {
                b = point;
                wb = point.residual;
                wa = if retained == -1 { 0.5 * wa } else { a.residual };
                retained = -1;
            }
        }
        let closest = if a.residual.abs() <= b.residual.abs() {
            a
        } else {
            b
        };
        Err(failure(self, closest, iterations))
    }
}
