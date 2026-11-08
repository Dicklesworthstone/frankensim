//! Implicit conduction between the original finite wall cells.
//!
//! An authored undirected link carries G (T_b - T_a) into a and the identical
//! negative heat into b. Here G = G_SI dt / C_f, with the SAME cell-capacity
//! scale C_f as FiniteThermalWall. Geometry, conductivity, contact resistance
//! and the correspondence between a link and a physical solid path remain
//! caller-owned. Links can resolve lateral conduction or declared connections
//! between wall planes; this does not create interior solid depth cells.
//!
//! After the actual scalar collision, the frozen fluid trace is g_out / w.
//! Eliminating the reflected pair gives K = 1 / (R + 1/(2w)). One backward-Euler
//! solve then advances every connected solid cell and its fluid exchange:
//! (C + K + L_G) delta_T = Q + K (T_trace - T_old) - L_G T_old.
//! Centered increments preserve constant equilibrium independently of datum.
//! Symmetric diagonal scaling precedes the existing sparse LDLT factorization.
//! The factor is prepared once; each step uses only retained linear workspace.
//!
//! This is fixed-capacity, fixed-resistance conduction on lumped wall cells,
//! with first-order time accuracy. It removes the explicit conduction CFL,
//! not the fluid solver's stability constraints or spatial/time refinement.
//! Wall values remain the sole physical temperature owner. All temperatures
//! and reflections here are speculative scratch until the caller publishes
//! each wall once. A failed solve or physical balance check publishes nothing.

use super::{FiniteThermalWall, FiniteThermalWallExchange, ThermalWall, WallReflection};
use fs_sparse::{Coo, DirectOrdering, LdltError, LdltFactor, LdltOptions, SymbolicLdlt};
use std::sync::Arc;

/// One existing normal link on either slab wall, in its original surface order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct WallCell {
    /// Zero is the lower wall; one is the upper wall.
    pub side: usize,
    /// x in the original D2Q5 slab; x + nx*y in the original D3Q7 slab.
    /// Native open carriers remove the selected normal coordinate and retain
    /// x-fast surface order.
    pub link: usize,
}

/// One undirected conductive path between existing finite solid cells.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallConductionLink {
    pub a: WallCell,
    pub b: WallCell,
    /// Finite nonnegative G_SI dt / C_f. Zero explicitly insulates.
    /// Parallel paths must be combined into one authored conductance.
    pub conductance: f64,
}

/// Preparation bounds and absolute local/plane first-law tolerance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallConductionControls {
    /// Maximum distinct endpoints of positive-conductance links.
    pub max_nodes: usize,
    /// Maximum authored links, including zero-conductance declarations.
    pub max_links: usize,
    /// Existing sparse numeric factor and scratch payload budget in bytes.
    pub max_factor_bytes: usize,
    /// Absolute heat tolerance in C_f * delta_T units for each original step.
    pub balance_tolerance: f64,
}

impl Default for WallConductionControls {
    fn default() -> Self {
        Self {
            max_nodes: 65_536,
            max_links: 262_144,
            max_factor_bytes: 256 * 1024 * 1024,
            balance_tolerance: 1e-10,
        }
    }
}

/// Admission and numerical refusals; accepted wall state is never repaired.
#[derive(Debug, Clone, PartialEq)]
pub enum WallConductionError {
    Invalid(&'static str),
    Capacity,
    /// A connected wall no longer has its prepared capacity and resistance.
    ChangedWall(WallCell),
    /// An authored state, trial endpoint, flux or storage is unrepresentable.
    Arithmetic(WallCell),
    /// None identifies the independently reduced whole-plane balance.
    Balance {
        cell: Option<WallCell>,
        defect: f64,
        tolerance: f64,
    },
    Sparse(LdltError),
}

impl core::fmt::Display for WallConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "wall conduction: {self:?}")
    }
}
impl std::error::Error for WallConductionError {}
impl From<LdltError> for WallConductionError {
    fn from(error: LdltError) -> Self {
        Self::Sparse(error)
    }
}

#[derive(Clone, Copy)]
struct Node {
    cell: WallCell,
    capacity: f64,
    resistance: f64,
    fluid_conductance: f64,
    scale: f64,
}

#[derive(Clone, Copy)]
struct Edge {
    a: usize,
    b: usize,
    conductance: f64,
}

struct Prepared {
    plane: usize,
    weight: f64,
    tolerance: f64,
    nodes: Vec<Node>,
    // Sorted unique cells make a complete plane a direct surface-index lookup.
    side_ranges: [[usize; 2]; 2],
    edges: Vec<Edge>,
    factor: Option<LdltFactor>,
}

// The same reciprocal edge value is reduced on its two endpoints. Compensation
// also keeps a large exchange from concealing small sources in the audit.
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, x: f64) {
        let next = self.value + x;
        self.correction += if self.value.abs() >= x.abs() {
            (self.value - next) + x
        } else {
            (x - next) + self.value
        };
        self.value = next;
    }
    fn total(self) -> f64 {
        self.value + self.correction
    }
}

#[derive(Clone, Copy, Default)]
struct Reflection {
    endpoint: f64,
    outgoing: f64,
    incoming: f64,
    surface: f64,
    audit: FiniteThermalWallExchange,
}

/// Prepared sparse solve and reusable scratch, with no independent solid state.
///
/// Cloning shares the immutable numeric factor. Scratch is copied so existing
/// slab snapshots can roll back together with populations and wall values.
#[derive(Clone)]
pub struct WallConduction {
    prepared: Arc<Prepared>,
    rhs: Vec<f64>,
    conducted: Vec<Sum>,
    reflections: Vec<Reflection>,
    ready: bool,
}

impl core::fmt::Debug for WallConduction {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("WallConduction")
            .field("plane", &self.plane())
            .field("normal_weight", &self.normal_weight())
            .field("nodes", &self.node_count())
            .field("links", &self.link_count())
            .finish_non_exhaustive()
    }
}

fn buffer<T: Clone + Default>(n: usize) -> Result<Vec<T>, WallConductionError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| WallConductionError::Capacity)?;
    values.resize(n, T::default());
    Ok(values)
}

fn solid(
    walls: &[Vec<ThermalWall>; 2],
    cell: WallCell,
) -> Result<FiniteThermalWall, WallConductionError> {
    let Some(ThermalWall::FiniteSolid(value)) =
        walls.get(cell.side).and_then(|side| side.get(cell.link))
    else {
        return Err(WallConductionError::Invalid(
            "a conduction endpoint must be an existing finite wall cell",
        ));
    };
    value
        .validate()
        .map_err(|_| WallConductionError::Arithmetic(cell))?;
    Ok(*value)
}

// A scaled difference can remain finite when subtracting opposite large
// temperature coordinates would overflow. No heat or temperature is clamped.
fn difference(coefficient: f64, first: f64, second: f64) -> f64 {
    let delta = first - second;
    if delta.is_finite() {
        coefficient * delta
    } else {
        coefficient * first - coefficient * second
    }
}

impl WallConduction {
    /// Prepare a fixed-property graph over existing wall cells.
    ///
    /// Only endpoints of positive links enter the factor; other links keep
    /// their original scalar wall update. Empty/zero-link graphs are valid.
    /// Repeated undirected pairs, self links, invalid endpoints and exceeded
    /// work bounds refuse. Changing a connected C or R requires preparation
    /// again; its temperature and prescribed heat may change each step.
    pub fn prepare(
        walls: &[Vec<ThermalWall>; 2],
        plane: usize,
        normal_weight: f64,
        links: &[WallConductionLink],
        controls: WallConductionControls,
    ) -> Result<Self, WallConductionError> {
        if plane == 0
            || walls.iter().any(|wall| wall.len() != plane)
            || !normal_weight.is_finite()
            || normal_weight <= 0.0
            || normal_weight > 0.5
            || controls.max_nodes == 0
            || controls.max_links == 0
            || controls.max_factor_bytes == 0
            || !controls.balance_tolerance.is_finite()
            || controls.balance_tolerance < 0.0
            || links.len() > controls.max_links
        {
            return Err(WallConductionError::Invalid(
                "invalid wall shape, weight or conduction controls",
            ));
        }
        let mut authored = Vec::new();
        authored
            .try_reserve_exact(links.len())
            .map_err(|_| WallConductionError::Capacity)?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(
                links
                    .len()
                    .checked_mul(2)
                    .ok_or(WallConductionError::Capacity)?,
            )
            .map_err(|_| WallConductionError::Capacity)?;
        for link in links {
            if link.a == link.b || !link.conductance.is_finite() || link.conductance < 0.0 {
                return Err(WallConductionError::Invalid(
                    "conduction requires distinct cells and finite nonnegative conductance",
                ));
            }
            solid(walls, link.a)?;
            solid(walls, link.b)?;
            let (a, b) = if link.a < link.b {
                (link.a, link.b)
            } else {
                (link.b, link.a)
            };
            authored.push(WallConductionLink {
                a,
                b,
                conductance: link.conductance,
            });
            if link.conductance > 0.0 {
                cells.extend_from_slice(&[a, b]);
            }
        }
        authored.sort_unstable_by_key(|link| (link.a, link.b));
        if authored
            .windows(2)
            .any(|pair| pair[0].a == pair[1].a && pair[0].b == pair[1].b)
        {
            return Err(WallConductionError::Invalid(
                "duplicate undirected wall conduction link",
            ));
        }
        cells.sort_unstable();
        cells.dedup();
        if cells.len() > controls.max_nodes {
            return Err(WallConductionError::Invalid(
                "wall conduction node budget exceeded",
            ));
        }
        let mut nodes = Vec::new();
        nodes
            .try_reserve_exact(cells.len())
            .map_err(|_| WallConductionError::Capacity)?;
        for &cell in &cells {
            let wall = solid(walls, cell)?;
            let fluid_conductance = 1.0 / (wall.resistance + 0.5 / normal_weight);
            if !fluid_conductance.is_finite() || fluid_conductance <= 0.0 {
                return Err(WallConductionError::Arithmetic(cell));
            }
            nodes.push(Node {
                cell,
                capacity: wall.heat_capacity,
                resistance: wall.resistance,
                fluid_conductance,
                scale: 0.0,
            });
        }
        let mut edges = Vec::new();
        edges
            .try_reserve_exact(authored.len())
            .map_err(|_| WallConductionError::Capacity)?;
        for link in authored.into_iter().filter(|link| link.conductance > 0.0) {
            edges.push(Edge {
                a: cells.binary_search(&link.a).expect("admitted endpoint"),
                b: cells.binary_search(&link.b).expect("admitted endpoint"),
                conductance: link.conductance,
            });
        }
        let mut diagonal = buffer::<Sum>(nodes.len())?;
        for (sum, node) in diagonal.iter_mut().zip(&nodes) {
            sum.add(node.capacity);
            sum.add(node.fluid_conductance);
        }
        for edge in &edges {
            diagonal[edge.a].add(edge.conductance);
            diagonal[edge.b].add(edge.conductance);
        }
        let mut matrix = Coo::new(nodes.len(), nodes.len());
        for (i, node) in nodes.iter_mut().enumerate() {
            let value = diagonal[i].total();
            node.scale = 1.0 / value.sqrt();
            if !value.is_finite() || value <= 0.0 || !node.scale.is_finite() || node.scale <= 0.0 {
                return Err(WallConductionError::Arithmetic(node.cell));
            }
            matrix.push(i, i, (value * node.scale) * node.scale);
        }
        for edge in &edges {
            let value = -(edge.conductance * nodes[edge.a].scale) * nodes[edge.b].scale;
            matrix.push(edge.a, edge.b, value);
            matrix.push(edge.b, edge.a, value);
        }
        let factor = if nodes.is_empty() {
            None
        } else {
            let matrix = matrix.assemble();
            let symbolic = SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd)?;
            let options = LdltOptions {
                max_factor_bytes: Some(controls.max_factor_bytes),
                ..LdltOptions::default()
            };
            let factor = symbolic.factor(&matrix, &options)?;
            if factor.inertia().positive != nodes.len() || factor.inertia().negative != 0 {
                return Err(WallConductionError::Invalid(
                    "wall conduction matrix lost positive definiteness",
                ));
            }
            Some(factor)
        };
        let n = nodes.len();
        let split = nodes.partition_point(|node| node.cell.side == 0);
        Ok(Self {
            prepared: Arc::new(Prepared {
                plane,
                weight: normal_weight,
                tolerance: controls.balance_tolerance,
                nodes,
                side_ranges: [[0, split], [split, n]],
                edges,
                factor,
            }),
            rhs: buffer(n)?,
            conducted: buffer(n)?,
            reflections: buffer(n)?,
            ready: false,
        })
    }

    #[must_use]
    pub fn plane(&self) -> usize {
        self.prepared.plane
    }
    #[must_use]
    pub fn normal_weight(&self) -> f64 {
        self.prepared.weight
    }
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.prepared.nodes.len()
    }
    #[must_use]
    pub fn link_count(&self) -> usize {
        self.prepared.edges.len()
    }

    /// Solve against the actual post-collision outgoing populations.
    ///
    /// No wall is mutated. A refusal invalidates all prior prepared
    /// reflections, even if the previous call succeeded. Each connected cell
    /// and each wall plane must close storage + fluid heat - source - internal
    /// conduction within the declared absolute tolerance before this returns.
    /// Factorization is reused and this method performs no heap allocation.
    pub(crate) fn prepare_reflections(
        &mut self,
        walls: &[Vec<ThermalWall>; 2],
        outgoing: impl Fn(usize, usize) -> f64,
    ) -> Result<(), WallConductionError> {
        self.ready = false;
        let prepared = &self.prepared;
        if walls.iter().any(|wall| wall.len() != prepared.plane) {
            return Err(WallConductionError::Invalid(
                "prepared wall conduction shape changed",
            ));
        }
        for node in &prepared.nodes {
            let wall = solid(walls, node.cell)?;
            if wall.heat_capacity != node.capacity || wall.resistance != node.resistance {
                return Err(WallConductionError::ChangedWall(node.cell));
            }
        }
        self.conducted.fill(Sum::default());
        for edge in &prepared.edges {
            let a = prepared.nodes[edge.a].cell;
            let b = prepared.nodes[edge.b].cell;
            let heat = difference(
                edge.conductance,
                solid(walls, b)?.temperature,
                solid(walls, a)?.temperature,
            );
            if !heat.is_finite() {
                return Err(WallConductionError::Arithmetic(a));
            }
            self.conducted[edge.a].add(heat);
            self.conducted[edge.b].add(-heat);
        }
        for (i, node) in prepared.nodes.iter().enumerate() {
            let wall = solid(walls, node.cell)?;
            let population = outgoing(node.cell.side, node.cell.link);
            let trace = population / prepared.weight;
            let mut rhs = self.conducted[i];
            rhs.add(wall.heat_per_step);
            rhs.add(difference(node.fluid_conductance, trace, wall.temperature));
            self.rhs[i] = node.scale * rhs.total();
            if !population.is_finite() || !trace.is_finite() || !self.rhs[i].is_finite() {
                return Err(WallConductionError::Arithmetic(node.cell));
            }
            self.reflections[i].outgoing = population;
        }
        if let Some(factor) = &prepared.factor {
            factor.solve_in_place(&mut self.rhs);
        }
        for (i, node) in prepared.nodes.iter().enumerate() {
            let old = solid(walls, node.cell)?.temperature;
            let endpoint = old + node.scale * self.rhs[i];
            if !endpoint.is_finite() {
                return Err(WallConductionError::Arithmetic(node.cell));
            }
            self.reflections[i].endpoint = endpoint;
        }
        self.conducted.fill(Sum::default());
        for edge in &prepared.edges {
            let heat = difference(
                edge.conductance,
                self.reflections[edge.b].endpoint,
                self.reflections[edge.a].endpoint,
            );
            if !heat.is_finite() {
                return Err(WallConductionError::Arithmetic(prepared.nodes[edge.a].cell));
            }
            self.conducted[edge.a].add(heat);
            self.conducted[edge.b].add(-heat);
        }
        let mut plane_accounts = [[Sum::default(); 4]; 2];
        for (i, node) in prepared.nodes.iter().enumerate() {
            let wall = solid(walls, node.cell)?;
            let reflection = &mut self.reflections[i];
            let requested_heat = difference(
                node.fluid_conductance,
                reflection.endpoint,
                reflection.outgoing / prepared.weight,
            );
            reflection.incoming = reflection.outgoing + requested_heat;
            reflection.surface =
                (0.5 * reflection.incoming + 0.5 * reflection.outgoing) / prepared.weight;
            let heat_in = reflection.incoming - reflection.outgoing;
            let stored_change = difference(node.capacity, reflection.endpoint, wall.temperature);
            let conducted_in = self.conducted[i].total();
            let source = wall.heat_per_step;
            let mut balance = Sum::default();
            for value in [stored_change, heat_in, -source, -conducted_in] {
                balance.add(value);
            }
            let defect = balance.total();
            if !reflection.incoming.is_finite()
                || !reflection.surface.is_finite()
                || !heat_in.is_finite()
                || !stored_change.is_finite()
                || !conducted_in.is_finite()
                || !defect.is_finite()
            {
                return Err(WallConductionError::Arithmetic(node.cell));
            }
            if defect.abs() > prepared.tolerance {
                return Err(WallConductionError::Balance {
                    cell: Some(node.cell),
                    defect,
                    tolerance: prepared.tolerance,
                });
            }
            reflection.audit = FiniteThermalWallExchange {
                heat_in,
                stored_change,
                source,
                conducted_in,
                maximum_link_defect: defect.abs(),
            };
            for (sum, value) in plane_accounts[node.cell.side].iter_mut().zip([
                stored_change,
                heat_in,
                -source,
                -conducted_in,
            ]) {
                sum.add(value);
            }
        }
        for account in plane_accounts {
            let mut balance = Sum::default();
            for sum in account {
                balance.add(sum.total());
            }
            let defect = balance.total();
            if !defect.is_finite() || defect.abs() > prepared.tolerance {
                return Err(WallConductionError::Balance {
                    cell: None,
                    defect,
                    tolerance: prepared.tolerance,
                });
            }
        }
        self.ready = true;
        Ok(())
    }

    /// Publish one accepted endpoint and reflection, once per original link.
    ///
    /// The slab owns this call sequence after successful preparation, with
    /// the same unmodified wall values and outgoing populations. Unconnected
    /// links use their existing scalar law, including independent finite cells.
    pub(crate) fn reflect(
        &self,
        side: usize,
        link: usize,
        wall: &mut ThermalWall,
        outgoing: f64,
        weight: f64,
    ) -> WallReflection {
        let [start, end] = self.prepared.side_ranges[side];
        let index = if end - start == self.prepared.plane {
            Some(start + link)
        } else {
            self.prepared.nodes[start..end]
                .binary_search_by_key(&link, |node| node.cell.link)
                .ok()
                .map(|offset| start + offset)
        };
        let Some(index) = index else {
            return wall.reflect_and_advance(outgoing, weight);
        };
        assert!(
            self.ready,
            "wall conduction must be prepared before streaming"
        );
        let reflection = self.reflections[index];
        debug_assert_eq!(weight, self.prepared.weight);
        debug_assert_eq!(outgoing, reflection.outgoing);
        let ThermalWall::FiniteSolid(solid) = wall else {
            unreachable!("prepared solid wall must remain finite until streaming");
        };
        solid.temperature = reflection.endpoint;
        WallReflection {
            incoming: reflection.incoming,
            temperature: reflection.surface,
            solid: Some(reflection.audit),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wall(temperature: f64, capacity: f64, resistance: f64, source: f64) -> ThermalWall {
        ThermalWall::FiniteSolid(FiniteThermalWall {
            temperature,
            heat_capacity: capacity,
            resistance,
            heat_per_step: source,
        })
    }
    fn temperature(wall: ThermalWall) -> f64 {
        let ThermalWall::FiniteSolid(value) = wall else {
            panic!("finite wall expected");
        };
        value.temperature
    }
    fn close(first: f64, second: f64, tolerance: f64) {
        assert!(
            (first - second).abs() <= tolerance,
            "{first:.17e} != {second:.17e}"
        );
    }
    fn link(a: (usize, usize), b: (usize, usize), conductance: f64) -> WallConductionLink {
        WallConductionLink {
            a: WallCell {
                side: a.0,
                link: a.1,
            },
            b: WallCell {
                side: b.0,
                link: b.1,
            },
            conductance,
        }
    }

    #[test]
    fn two_cell_coupling_matches_analytic_storage_and_reciprocal_heat() {
        for weight in [1.0 / 6.0, 0.125] {
            let mut walls = [
                vec![wall(0.9, 1.7, 0.6, 0.07)],
                vec![wall(0.2, 4.1, 3.4, -0.02)],
            ];
            let original = walls.clone();
            let conductance = 1.8;
            let mut graph = WallConduction::prepare(
                &walls,
                1,
                weight,
                &[link((0, 0), (1, 0), conductance)],
                WallConductionControls::default(),
            )
            .unwrap();
            let fluid = [0.3, 0.6];
            let k = [1.0 / (0.6 + 0.5 / weight), 1.0 / (3.4 + 0.5 / weight)];
            let diagonal = [1.7 + k[0] + conductance, 4.1 + k[1] + conductance];
            let rhs = [
                1.7 * 0.9 + 0.07 + k[0] * fluid[0],
                4.1 * 0.2 - 0.02 + k[1] * fluid[1],
            ];
            let determinant = diagonal[0] * diagonal[1] - conductance * conductance;
            let expected = [
                (rhs[0] * diagonal[1] + conductance * rhs[1]) / determinant,
                (rhs[1] * diagonal[0] + conductance * rhs[0]) / determinant,
            ];
            graph
                .prepare_reflections(&walls, |side, _| weight * fluid[side])
                .unwrap();
            assert_eq!(walls, original, "preparation must not publish state");
            let mut audits = [FiniteThermalWallExchange::default(); 2];
            for side in 0..2 {
                let reflection =
                    graph.reflect(side, 0, &mut walls[side][0], weight * fluid[side], weight);
                audits[side] = reflection.solid.unwrap();
                close(temperature(walls[side][0]), expected[side], 3e-15);
                let resistance = if side == 0 { 0.6 } else { 3.4 };
                close(
                    expected[side] - reflection.temperature,
                    resistance * audits[side].heat_in,
                    3e-15,
                );
                close(audits[side].defect(), 0.0, 3e-15);
            }
            assert!(audits[0].conducted_in < 0.0);
            assert_eq!(audits[0].conducted_in, -audits[1].conducted_in);
            close(
                audits
                    .iter()
                    .map(|a| a.stored_change + a.heat_in - a.source)
                    .sum(),
                0.0,
                3e-15,
            );
        }
    }

    #[test]
    fn lateral_fourier_mode_converges_in_space_and_stiff_conduction_stays_bounded() {
        for weight in [1.0 / 6.0, 0.125] {
            let mut errors = Vec::new();
            for count in [8usize, 16, 32] {
                let dx = 1.0 / count as f64;
                let dt = 0.03;
                let resistance = 1e12;
                let mut walls = [
                    (0..count)
                        .map(|i| {
                            wall(
                                1.0 + 0.2 * (std::f64::consts::TAU * i as f64 * dx).cos(),
                                dx,
                                resistance,
                                0.0,
                            )
                        })
                        .collect::<Vec<_>>(),
                    vec![ThermalWall::Adiabatic; count],
                ];
                let links = (0..count)
                    .map(|i| link((0, i), (0, (i + 1) % count), dt / dx))
                    .collect::<Vec<_>>();
                let mut graph = WallConduction::prepare(
                    &walls,
                    count,
                    weight,
                    &links,
                    WallConductionControls::default(),
                )
                .unwrap();
                graph.prepare_reflections(&walls, |_, _| weight).unwrap();
                let mut ledger = Sum::default();
                let mut conducted = Sum::default();
                let mut maximum_error = 0.0f64;
                let continuum_amplitude = 0.2 / (1.0 + dt * std::f64::consts::TAU.powi(2));
                for i in 0..count {
                    let response = graph.reflect(0, i, &mut walls[0][i], weight, weight);
                    let audit = response.solid.unwrap();
                    ledger.add(audit.stored_change + audit.heat_in);
                    conducted.add(audit.conducted_in);
                    let t = temperature(walls[0][i]);
                    assert!((0.8..=1.2).contains(&t));
                    let expected =
                        1.0 + continuum_amplitude * (std::f64::consts::TAU * i as f64 * dx).cos();
                    maximum_error = maximum_error.max((t - expected).abs());
                }
                close(ledger.total(), 0.0, 2e-14);
                close(conducted.total(), 0.0, 2e-14);
                errors.push(maximum_error);
            }
            assert!(
                errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
                "{errors:?}"
            );
            let mut walls = [
                vec![wall(1.0, 1.0, 0.0, 0.0), wall(0.0, 1.0, 0.0, 0.0)],
                vec![ThermalWall::Adiabatic; 2],
            ];
            let mut graph = WallConduction::prepare(
                &walls,
                2,
                weight,
                &[link((0, 0), (0, 1), 1_000.0)],
                WallConductionControls::default(),
            )
            .unwrap();
            graph
                .prepare_reflections(&walls, |_, _| 0.5 * weight)
                .unwrap();
            for i in 0..2 {
                let response = graph.reflect(0, i, &mut walls[0][i], 0.5 * weight, weight);
                assert!(response.solid.unwrap().maximum_link_defect < 1e-10);
                assert!((0.0..=1.0).contains(&temperature(walls[0][i])));
            }
            assert!((temperature(walls[0][0]) - temperature(walls[0][1])).abs() < 0.001);
        }
    }

    #[test]
    fn equilibrium_scaling_and_unconnected_cells_preserve_original_storage_owner() {
        for weight in [1.0 / 6.0, 0.125] {
            let mut walls = [
                vec![
                    wall(0.75, 1.0, 0.0, 0.0),
                    wall(0.75, 1e20, 0.0, 0.0),
                    wall(0.4, 2.0, 3.0, 0.07),
                ],
                vec![ThermalWall::Adiabatic; 3],
            ];
            let original = walls.clone();
            let mut graph = WallConduction::prepare(
                &walls,
                3,
                weight,
                &[link((0, 0), (0, 1), 1.0)],
                WallConductionControls::default(),
            )
            .unwrap();
            graph
                .prepare_reflections(&walls, |_, _| 0.75 * weight)
                .unwrap();
            for i in 0..2 {
                let response = graph.reflect(0, i, &mut walls[0][i], 0.75 * weight, weight);
                assert_eq!(
                    response.solid.unwrap(),
                    FiniteThermalWallExchange::default()
                );
                assert_eq!(walls[0][i], original[0][i]);
            }
            let mut expected = original[0][2];
            let scalar = expected.reflect_and_advance(0.75 * weight, weight);
            let response = graph.reflect(0, 2, &mut walls[0][2], 0.75 * weight, weight);
            assert_eq!(walls[0][2], expected);
            assert_eq!(response.incoming, scalar.incoming);
            assert_eq!(response.solid, scalar.solid);
            let mut empty =
                WallConduction::prepare(&walls, 3, weight, &[], WallConductionControls::default())
                    .unwrap();
            assert_eq!(empty.node_count(), 0);
            empty
                .prepare_reflections(&walls, |_, _| 0.75 * weight)
                .unwrap();
        }
    }

    #[test]
    fn failed_balance_or_changed_properties_never_publish_candidate() {
        for weight in [1.0 / 6.0, 0.125] {
            let mut walls = [
                vec![wall(1.0, 1e300, 0.0, 0.0), wall(1.0, 1.0, 0.0, 0.0)],
                vec![ThermalWall::Adiabatic; 2],
            ];
            let mut graph = WallConduction::prepare(
                &walls,
                2,
                weight,
                &[link((0, 0), (0, 1), 1.0)],
                WallConductionControls::default(),
            )
            .unwrap();
            let original = walls.clone();
            graph.prepare_reflections(&walls, |_, _| weight).unwrap();
            assert!(graph.ready);
            assert!(matches!(
                graph.prepare_reflections(&walls, |_, _| 0.25 * weight),
                Err(WallConductionError::Balance { .. })
            ));
            assert!(!graph.ready);
            assert_eq!(walls, original);
            let ThermalWall::FiniteSolid(solid) = &mut walls[0][1] else {
                unreachable!();
            };
            solid.resistance = 2.0;
            let changed = walls.clone();
            assert!(matches!(
                graph.prepare_reflections(&walls, |_, _| weight),
                Err(WallConductionError::ChangedWall(_))
            ));
            assert_eq!(walls, changed);
        }
    }

    #[test]
    fn invalid_graph_and_numeric_factor_budget_refuse_before_publication() {
        let walls = [
            vec![wall(1.0, 1.0, 0.0, 0.0), wall(0.0, 1.0, 0.0, 0.0)],
            vec![ThermalWall::Adiabatic; 2],
        ];
        let controls = WallConductionControls::default();
        for links in [
            vec![link((0, 0), (0, 0), 1.0)],
            vec![link((0, 0), (1, 0), 1.0)],
            vec![link((0, 0), (0, 2), 1.0)],
            vec![link((0, 0), (0, 1), f64::NAN)],
            vec![link((0, 0), (0, 1), 1.0), link((0, 1), (0, 0), 2.0)],
        ] {
            assert!(WallConduction::prepare(&walls, 2, 0.125, &links, controls).is_err());
        }
        assert!(
            WallConduction::prepare(
                &walls,
                2,
                0.125,
                &[link((0, 0), (0, 1), 1.0)],
                WallConductionControls {
                    max_nodes: 1,
                    ..controls
                }
            )
            .is_err()
        );
        assert!(matches!(
            WallConduction::prepare(
                &walls,
                2,
                0.125,
                &[link((0, 0), (0, 1), 1.0)],
                WallConductionControls {
                    max_factor_bytes: 1,
                    ..controls
                }
            ),
            Err(WallConductionError::Sparse(_))
        ));
    }
}
