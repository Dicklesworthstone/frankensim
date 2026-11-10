//! Exact frozen affine current offsets on the original circuit descriptor.
//! These internal ports lower I=G V+b; they are not independent heat sources.
use super::*;

pub struct NortonPorts {
    pub graph: CircuitGraph,
    pub values: Vec<f64>,
    authored: usize,
}

pub struct NortonWork {
    pub authored_j: f64,
    pub offset_j: f64,
    pub total_j: f64,
}

impl NortonPorts {
    /// Offsets are signed INTO the conductor in original terminal order.
    /// Each non-reference offset is a source from node i to node zero; the
    /// reference offset follows conserved current, without a grounding choice.
    pub fn new(
        graph: &CircuitGraph,
        nodes: &[usize],
        offset_current_a: &[f64],
        ports: &[f64],
        current_tolerance_a: f64,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error> {
        check(cancel)?;
        if nodes.len() < 2
            || nodes.len() != offset_current_a.len()
            || nodes.len() > 33
            || nodes.iter().any(|&node| node >= graph.node_count)
            || nodes
                .iter()
                .enumerate()
                .any(|(i, node)| nodes[..i].contains(node))
            || !current_tolerance_a.is_finite()
            || current_tolerance_a <= 0.0
            || ports
                .iter()
                .chain(offset_current_a)
                .any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite affine terminal currents and source ports required",
            ));
        }
        let mut authored = 0usize;
        for &(_, _, branch) in &graph.branches {
            check(cancel)?;
            if let Branch::VoltageSource { port } | Branch::CurrentSource { port } = branch {
                authored = authored.max(
                    port.checked_add(1)
                        .ok_or_else(|| invalid("authored circuit port count overflows"))?,
                );
            }
        }
        if authored != ports.len()
            || authored
                .checked_add(nodes.len() - 1)
                .is_none_or(|count| count > 32)
            || sum(offset_current_a.iter().copied())?.abs() > current_tolerance_a
        {
            return Err(invalid(
                "bounded complete authored ports and conserved motional currents required",
            ));
        }
        let mut graph = graph.clone();
        let mut values = ports.to_vec();
        for (&node, &offset) in nodes[1..].iter().zip(&offset_current_a[1..]) {
            check(cancel)?;
            graph
                .branches
                .push((nodes[0], node, Branch::CurrentSource { port: values.len() }));
            values.push(offset);
        }
        check(cancel)?;
        Ok(Self {
            graph,
            values,
            authored,
        })
    }

    /// Independently pair original authored ports and the exact affine stamp.
    /// Physical Joule/Lorentz work is recovered from the native field, never
    /// inferred from this internal port's work or homogeneous circuit loss.
    pub fn work(&self, outputs: &[f64], duration_s: f64) -> Result<NortonWork, Error> {
        if outputs.len() != self.values.len()
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || outputs.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                "complete finite circuit outputs and positive interval required",
            ));
        }
        let products = || self.values.iter().zip(outputs).map(|(&u, &y)| u * y);
        let authored_j = finite(duration_s * sum(products().take(self.authored))?)?;
        let offset_j = finite(duration_s * sum(products().skip(self.authored))?)?;
        let total_j = finite(duration_s * sum(products())?)?;
        Ok(NortonWork {
            authored_j,
            offset_j,
            total_j,
        })
    }
}

fn sum(values: impl IntoIterator<Item = f64>) -> Result<f64, Error> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        finite(value)?;
        let next = finite(total + value)?;
        correction = finite(
            correction
                + if total.abs() >= value.abs() {
                    (total - next) + value
                } else {
                    (value - next) + total
                },
        )?;
        total = next;
    }
    Ok(finite(total + correction)?)
}
