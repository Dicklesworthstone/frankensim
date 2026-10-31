//! Steady conservative finite-volume energy equation over fluid and solid
//! voxels.
//!
//! Every face contributes one outward total energy flux written as
//! `J = d T_P - a T_N - r` (Patankar, *Numerical Heat Transfer and Fluid
//! Flow*, 1980, §5.2–5.4): for an interior face with outward heat-capacity
//! flux `F = rho c_p Q` and diffusive conductance `D = A k_h / dx` (on a
//! graded grid `dx` is the centre distance and `k_h` weights the half widths),
//! `a = D A(|F|/D) + max(-F, 0)` and `d = a + F`; the two sides of one face
//! therefore carry equal and opposite `J` for ANY temperatures, so the global
//! energy balance telescopes exactly. `k_h = 2 k_P k_N / (k_P + k_N)` is the
//! series (harmonic) conductance of two half cells, exact for
//! piecewise-constant conductivity meeting at the face.
//! With authored k(T) tables, each bulk half cell integrates its own law in
//! temperature. One common flux satisfies both surface temperatures and the
//! declared contact jump. Positive secants preserve the conservative stencil;
//! a Robin surface solves the integrated half cell in series with its film.

use fs_exec::CancelGate;
use fs_sparse::Coo;

use super::conductance::{ConductivityCurve, interface_conductance};
use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::flow::{FlowField, scale_rows};
use super::krylov::bicgstab_ilu0;
use super::{ChtError, finite, finite_positive};
use crate::d3q19::Face3;

mod material_field;
mod enthalpy_field;
pub use enthalpy_field::{AdvectiveEnthalpyField, AffineSpecificEnthalpy};
mod face_contact;
pub use face_contact::FaceContactResistance;

/// Thermal rule on one domain face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalFace {
    /// Zero heat flux. Refuses if flow crosses the face.
    Adiabatic,
    /// Prescribed temperature (K) on the face plane, half a cell from the
    /// boundary cell centre.
    Temperature(f64),
    /// Heat flux INTO the domain, W/m^2.
    HeatFlux(f64),
    /// Film coefficient `h` (W/(m^2 K)) to an ambient temperature (K), in
    /// series with the boundary half cell.
    Convective {
        /// Film coefficient, W/(m^2 K).
        h: f64,
        /// Ambient temperature, K.
        ambient: f64,
    },
    /// Open inflow face: incoming fluid carries `temperature`; the face is
    /// also a diffusive Dirichlet boundary at that temperature. Solid cells
    /// on this face are treated as adiabatic.
    Inflow {
        /// Inflow temperature, K.
        temperature: f64,
    },
    /// Open outflow face: outgoing fluid carries its cell temperature, with
    /// zero diffusive flux. Any backflow carries `backflow_temperature`.
    /// Solid cells on this face are treated as adiabatic.
    Outflow {
        /// Temperature carried by any reversed flow, K.
        backflow_temperature: f64,
    },
}

impl ThermalFace {
    const fn is_open(self) -> bool {
        matches!(self, Self::Inflow { .. } | Self::Outflow { .. })
    }
}

/// Thermal boundary conditions, heat sources, and fixed-temperature cells.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalSetup {
    /// One rule per domain face in [`Face3::ALL`] order.
    pub faces: [ThermalFace; 6],
    /// Heat generated in each cell, W. Empty means no sources; otherwise one
    /// entry per cell.
    pub power_w: Vec<f64>,
    /// Cells held at a prescribed temperature (cell index, K). Their implied
    /// heat injection is reported in [`EnergyBalance::fixed_cell_injection_w`].
    pub fixed_temperature: Vec<(usize, f64)>,
    /// Thermal contact (interface) resistances between pairs of solid
    /// materials, in series on every face the two materials share.
    pub contacts: Vec<ContactResistance>,
    /// Explicit face-local resistances, including joints between cells of the
    /// SAME material. A face cannot also have a material-pair contact rule.
    /// These interfaces have no storage and admit no through-flow; they do not
    /// change either cell's mechanical occupancy or material assignment.
    pub face_contacts: Vec<FaceContactResistance>,
    /// Linear exchanges `G (T_cell - T_sink)` leaving individual cells (a
    /// linearized radiation exchange, a compact thermal model's link).
    pub cell_sinks: Vec<CellSink>,
    /// Turbulent conductivity added to each FLUID cell, W/(m K) (empty: none;
    /// otherwise one entry per cell, for example `FvFlow::eddy_conductivity`).
    pub eddy_conductivity_w_m_k: Vec<f64>,
    /// Two-resistor compact thermal models of components.
    pub compact_components: Vec<CompactComponent>,
}

/// A JEDEC two-resistor compact thermal model (JESD15-3) of a packaged
/// component occupying the solid cells `lo..hi`: the dissipated power enters
/// one junction node, which reaches the case top through
/// `junction_to_case_k_w` and the board through `junction_to_board_k_w`,
/// each resistor spread over its face of the box in proportion to area; the
/// sides are adiabatic and the box interior is collapsed (its cells report
/// the junction temperature and carry no conduction of their own). The box
/// blocks flow like any solid. Both faces must border cells inside the
/// domain that belong to no other component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CompactComponent {
    /// First cells of the box (inclusive).
    pub lo: [usize; 3],
    /// Last cells of the box (exclusive).
    pub hi: [usize; 3],
    /// The face of the box that sits on the board; its opposite is the case
    /// top.
    pub board_face: Face3,
    /// Dissipated power, W.
    pub power_w: f64,
    /// Junction-to-case (top) resistance, K/W.
    pub junction_to_case_k_w: f64,
    /// Junction-to-board resistance, K/W.
    pub junction_to_board_k_w: f64,
}

impl CompactComponent {
    /// Area of the case-top (= board) face, m^2.
    fn face_area(&self, domain: &VoxelDomain) -> f64 {
        let axis = self.board_face as usize / 2;
        (0..3)
            .filter(|&a| a != axis)
            .map(|a| {
                (self.lo[a]..self.hi[a])
                    .map(|i| domain.width(a, i))
                    .sum::<f64>()
            })
            .product()
    }
}

/// Steady state of one compact component.
#[derive(Debug, Clone, PartialEq)]
pub struct JunctionSolution {
    /// Junction temperature, K.
    pub temperature_k: f64,
    /// Heat leaving through the junction-to-case resistor, W.
    pub case_w: f64,
    /// Heat leaving through the junction-to-board resistor, W.
    pub board_w: f64,
    /// The component's cell box (`lo..hi`).
    pub cells: ([usize; 3], [usize; 3]),
    /// Outside cells the junction couples to, with their conductances, W/K
    /// (the resistor's area share in series with the cell's half width).
    pub links: Vec<(usize, f64)>,
}

/// A linear heat path from one cell to a fixed temperature: the heat
/// leaving the cell is `conductance_w_k * (T_cell - temperature_k)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellSink {
    /// Cell index.
    pub cell: usize,
    /// Conductance, W/K (non-negative).
    pub conductance_w_k: f64,
    /// Sink temperature, K.
    pub temperature_k: f64,
}

/// Per-area thermal resistance of the interface between two solid
/// materials (a thermal interface material, a bonded joint, a pressed
/// contact), `R'' = dT / q''` in m^2 K / W.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactResistance {
    /// The two material indices into the `solids` slice (distinct).
    pub materials: (u16, u16),
    /// Per-area resistance, m^2 K / W (non-negative).
    pub resistance_m2_k_w: f64,
}

impl ThermalSetup {
    /// Face rules only; no sources or fixed cells.
    #[must_use]
    pub fn new(faces: [ThermalFace; 6]) -> Self {
        Self {
            faces,
            power_w: Vec::new(),
            fixed_temperature: Vec::new(),
            contacts: Vec::new(),
            face_contacts: Vec::new(),
            cell_sinks: Vec::new(),
            eddy_conductivity_w_m_k: Vec::new(),
            compact_components: Vec::new(),
        }
    }

    /// Spread `total_w` uniformly (by volume) over the cells whose centre
    /// satisfies `inside`; returns the number of heated cells.
    ///
    /// # Panics
    /// If `power_w` is non-empty with the wrong length.
    pub fn add_uniform_power(
        &mut self,
        domain: &VoxelDomain,
        total_w: f64,
        mut inside: impl FnMut([f64; 3]) -> bool,
    ) -> usize {
        if self.power_w.is_empty() {
            self.power_w = vec![0.0; domain.cell_count()];
        }
        assert_eq!(
            self.power_w.len(),
            domain.cell_count(),
            "power_w length mismatch"
        );
        let cells: Vec<usize> = (0..domain.cell_count())
            .filter(|&c| {
                let [x, y, z] = domain.coords(c);
                inside(domain.center(x, y, z))
            })
            .collect();
        // By volume: uniform power density on a graded grid as well.
        let volume: f64 = cells.iter().map(|&c| domain.volume(c)).sum();
        if volume > 0.0 {
            for &c in &cells {
                self.power_w[c] += total_w * domain.volume(c) / volume;
            }
        }
        cells.len()
    }
}

/// Convection discretization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConvectionScheme {
    /// Patankar's power-law scheme: `A(P) = max(0, (1 - 0.1 |P|)^5)`;
    /// close to the exact 1-D exponential profile at any cell Péclet number.
    #[default]
    PowerLaw,
    /// First-order upwind: `A(P) = 1`.
    Upwind,
}

/// Energy solve configuration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyConfig {
    /// Convection scheme.
    pub scheme: ConvectionScheme,
    /// Relative residual of the Jacobi-row-scaled system.
    pub tolerance: f64,
    /// Krylov iteration budget.
    pub max_iterations: usize,
}

impl Default for EnergyConfig {
    fn default() -> Self {
        Self {
            scheme: ConvectionScheme::PowerLaw,
            tolerance: 1e-12,
            max_iterations: 50_000,
        }
    }
}

/// Global energy accounting, W. Positive `boundary_outflow_w` leaves the
/// domain. `residual_w = source_w + fixed_cell_injection_w -
/// boundary_outflow_w` is computed from the boundary and fixed-cell fluxes
/// alone, independently of the solver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyBalance {
    /// Heat generated in free cells.
    pub source_w: f64,
    /// Heat the fixed-temperature cells inject to hold their temperature.
    pub fixed_cell_injection_w: f64,
    /// Net total energy flux (advective + diffusive) out of the domain.
    pub boundary_outflow_w: f64,
    /// Net upwind advective part of `boundary_outflow_w` (informational;
    /// constant-Cp mode references enthalpy to 0 K; an explicit caloric field
    /// uses its common supplied datum. Net boundary mass flow must vanish for
    /// this flux to be independent of that datum).
    pub advective_outflow_w: f64,
    /// Heat leaving through cell sinks (`ThermalSetup::cell_sinks`).
    pub sink_outflow_w: f64,
    /// `source + fixed - outflow - sink`.
    pub residual_w: f64,
    /// `|residual| / (|source| + |fixed| + sum |boundary flux|)`.
    pub relative_residual: f64,
}

/// Solver and discretization evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct EnergyReport {
    /// Unknowns (every cell).
    pub unknowns: usize,
    /// Stored nonzeros.
    pub nonzeros: usize,
    /// BiCGStab iterations.
    pub iterations: usize,
    /// Recomputed relative residual of the row-scaled system.
    pub relative_residual: f64,
    /// Largest interior-face cell Péclet number `|F| / D`.
    pub max_cell_peclet: f64,
    /// Convection scheme used.
    pub scheme: ConvectionScheme,
    /// Global energy accounting.
    pub balance: EnergyBalance,
}

/// Temperature field and its evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct EnergySolution {
    /// Temperature per cell, K.
    pub temperature: Vec<f64>,
    /// Cell-centre conductivity along x, y, z, W/(m K), retained for local
    /// material reporting. Authored k(T) tables integrate the bulk half cells;
    /// these local values alone do not reconstruct nonlinear face coefficients.
    pub conductivity: Vec<[f64; 3]>,
    /// Canonical `(lower cell, upper cell)` face resistances actually used by
    /// this solve, m² K/W. Kept for consistent downstream heat-flow observers.
    pub face_contact_resistance_m2_k_w: std::collections::BTreeMap<(usize, usize), f64>,
    /// Ordinary solid/fluid faces as `(solid cell, fluid cell, conductance W/K)`.
    /// These are the actual accepted face coefficients, including nonlinear
    /// bulk response and contact. Compact junctions retain their own links.
    /// Storage grows with boundary faces, not with all grid faces.
    pub solid_fluid_links: Vec<(usize, usize, f64)>,
    /// One entry per `ThermalSetup::compact_components`, in order.
    pub junctions: Vec<JunctionSolution>,
    /// Solve evidence.
    pub report: EnergyReport,
}

impl EnergySolution {
    /// Largest temperature over cells satisfying `select`, with its cell.
    #[must_use]
    pub fn max_where(&self, mut select: impl FnMut(usize) -> bool) -> Option<(usize, f64)> {
        let mut best: Option<(usize, f64)> = None;
        for (c, &t) in self.temperature.iter().enumerate() {
            if select(c) && best.is_none_or(|(_, b)| t > b) {
                best = Some((c, t));
            }
        }
        best
    }

    /// Flux-weighted (mixed-mean) fluid temperature of the cell layer `x`,
    /// using the mean of the layer's two `x`-face fluxes as the weight.
    /// `None` when the layer carries no net streamwise flow.
    #[must_use]
    pub fn bulk_temperature_x(
        &self,
        domain: &VoxelDomain,
        flow: &FlowField,
        x: usize,
    ) -> Option<f64> {
        let [_, ny, nz] = domain.dims();
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for z in 0..nz {
            for y in 0..ny {
                let c = domain.index(x, y, z);
                if !domain.is_fluid(c) {
                    continue;
                }
                let w = 0.5 * (flow.flux_x(x, y, z) + flow.flux_x(x + 1, y, z));
                num = w.mul_add(self.temperature[c], num);
                den += w;
            }
        }
        (den != 0.0).then(|| num / den)
    }

    /// Heat conducted from solid cells into fluid cells across every
    /// fluid/solid face, W, including any explicitly resolved face resistance.
    #[must_use]
    pub fn solid_to_fluid_heat_w(&self, domain: &VoxelDomain) -> f64 {
        let mut total = self.solid_fluid_links.iter().map(|&(solid, fluid, g)| {
            g * (self.temperature[solid] - self.temperature[fluid])
        }).sum();
        for junction in &self.junctions {
            for &(n, g) in &junction.links {
                if domain.is_fluid(n) {
                    total += g * (junction.temperature_k - self.temperature[n]);
                }
            }
        }
        total
    }
}

/// One face's contribution `J = diag T_P - off T_N - rhs`.
#[derive(Debug, Clone, Copy)]
struct FaceTerm {
    diag: f64,
    off: Option<(usize, f64)>,
    rhs: f64,
    /// Upwind advective part for the balance report (boundary faces only).
    advective: Advective,
}

#[derive(Debug, Clone, Copy)]
enum Advective {
    None,
    OwnCell(f64),
    Affine { mass_flow_kg_s: f64, enthalpy: AffineSpecificEnthalpy },
    Fixed(f64),
}

struct Context<'a> {
    domain: &'a VoxelDomain,
    flow: &'a FlowField,
    solids: &'a [SolidMaterial],
    /// Previous Picard state, only for the authored k(T) table path.
    conductivity_temperature: Option<&'a [f64]>,
    faces: [ThermalFace; 6],
    /// Conductivity along x, y, z per cell.
    k: Vec<[f64; 3]>,
    /// Contact resistances keyed by the ordered material pair.
    contacts: Vec<((u16, u16), f64)>,
    face_contacts: std::collections::BTreeMap<(usize, usize), f64>,
    rho_c: f64,
    density: f64,
    enthalpy: Option<AdvectiveEnthalpyField<'a>>,
    scheme: ConvectionScheme,
}

impl Context<'_> {
    fn checked_term(&self, c: usize, f: usize) -> Result<FaceTerm, ChtError> {
        let term = self.term(c, f);
        if !term.diag.is_finite()
            || !term.rhs.is_finite()
            || term
                .off
                .is_some_and(|(_, coefficient)| !coefficient.is_finite())
        {
            return Err(ChtError::InvalidInput {
                field: "thermal.face_coefficients",
                reason: "face energy coefficients or nonlinear interface are unrepresentable".into(),
            });
        }
        Ok(term)
    }

    fn contact(&self, c: usize, n: usize) -> f64 {
        if let Some(&resistance) = self.face_contacts.get(&(c.min(n), c.max(n))) {
            return resistance;
        }
        if self.contacts.is_empty() {
            return 0.0;
        }
        match (self.domain.voxel_at(c), self.domain.voxel_at(n)) {
            (Voxel::Solid(a), Voxel::Solid(b)) if a != b => {
                let key = (a.min(b), a.max(b));
                self.contacts
                    .iter()
                    .find(|(pair, _)| *pair == key)
                    .map_or(0.0, |(_, r)| *r)
            }
            _ => 0.0,
        }
    }

    fn conductivity_curve(&self, c: usize, axis: usize) -> ConductivityCurve<'_> {
        if self.conductivity_temperature.is_some()
            && let Voxel::Solid(m) = self.domain.voxel_at(c)
        {
            let material = &self.solids[usize::from(m)];
            if !material.conductivity_table.is_empty() {
                return ConductivityCurve::tabulated(
                    material.conductivity_w_m_k,
                    material.axis_conductivity()[axis],
                    &material.conductivity_table,
                );
            }
        }
        // Resolved molecular fields and fluid eddy conductivity are already
        // frozen by their owners. No temperature curve is inferred from them.
        ConductivityCurve::constant(self.k[c][axis])
    }

    fn diffusive_conductance(&self, c: usize, n: usize, axis: usize) -> f64 {
        let area = self.domain.face_area(c, axis);
        let half_c = 0.5 * self.domain.widths(c)[axis];
        let half_n = 0.5 * self.domain.widths(n)[axis];
        let contact = self.contact(c, n);
        let has_curve = |cell| match self.domain.voxel_at(cell) {
            Voxel::Solid(m) => !self.solids[usize::from(m)].conductivity_table.is_empty(),
            Voxel::Fluid => false,
        };
        let d = if let Some(t) = self.conductivity_temperature
            && (has_curve(c) || has_curve(n))
        {
            if let (Voxel::Solid(a), Voxel::Solid(b)) =
                (self.domain.voxel_at(c), self.domain.voxel_at(n))
                && a == b && contact == 0.0
            {
                // One continuous curve has a direct Kirchhoff secant.
                let solid = &self.solids[usize::from(a)];
                area * solid.axis_conductivity()[axis]
                    * solid.conductivity_secant_factor(t[c], t[n])
                    / (half_c + half_n)
            } else {
                // Both rows use exactly the same ordering and positive
                // conductance. Each bulk has its own integrated curve; the
                // two surface temperatures retain the declared contact jump.
                let (a, b) = (c.min(n), c.max(n));
                area * interface_conductance(
                    self.conductivity_curve(a, axis), t[a],
                    0.5 * self.domain.widths(a)[axis],
                    self.conductivity_curve(b, axis), t[b],
                    0.5 * self.domain.widths(b)[axis], contact,
                ).unwrap_or(f64::NAN)
            }
        } else {
            area / (half_c / self.k[c][axis] + contact + half_n / self.k[n][axis])
        };
        // A failed or unrepresentable face must not silently become an
        // insulator. checked_term gives it the existing typed refusal.
        if d.is_finite() && d > 0.0 { d } else { f64::NAN }
    }

    fn weight(&self, peclet: f64) -> f64 {
        match self.scheme {
            ConvectionScheme::Upwind => 1.0,
            ConvectionScheme::PowerLaw => {
                let t = 0.1f64.mul_add(-peclet.abs(), 1.0).max(0.0);
                let t2 = t * t;
                t2 * t2 * t
            }
        }
    }

    fn enthalpy_boundary(&self, c: usize, f: usize) -> FaceTerm {
        let field = self.enthalpy.expect("enthalpy path");
        let mass = self.density * self.flow.outward(self.domain, c, f);
        if mass >= 0.0 {
            let enthalpy = field.cells[c].expect("admitted fluid boundary");
            FaceTerm {
                diag: mass * enthalpy.heat_capacity_j_kg_k, off: None,
                rhs: -mass * enthalpy.intercept(),
                advective: Advective::Affine { mass_flow_kg_s: mass, enthalpy },
            }
        } else {
            let j = mass * field.boundary_j_kg[f].expect("admitted boundary enthalpy");
            FaceTerm { diag: 0.0, off: None, rhs: -j, advective: Advective::Fixed(j) }
        }
    }

    fn term(&self, c: usize, f: usize) -> FaceTerm {
        let axis = f / 2;
        let area = self.domain.face_area(c, axis);
        let half_c = 0.5 * self.domain.widths(c)[axis];
        let flux = self.rho_c * self.flow.outward(self.domain, c, f);
        let kc = self.k[c][f / 2];
        if let Some(n) = self.domain.neighbor(c, f) {
            let d = self.diffusive_conductance(c, n, axis);
            if let Some(field) = self.enthalpy {
                let mass = self.density * self.flow.outward(self.domain, c, f);
                // Solid faces are impermeable. At a fluid face both rows use
                // exactly the same donor tangent, including its affine offset.
                if mass != 0.0 {
                    let upstream = if mass > 0.0 { c } else { n };
                    let h = field.cells[upstream].expect("admitted fluid donor");
                    let rate = mass * h.heat_capacity_j_kg_k;
                    return FaceTerm {
                        diag: d + rate.max(0.0),
                        off: Some((n, d + (-rate).max(0.0))),
                        rhs: -mass * h.intercept(),
                        advective: Advective::None,
                    };
                }
            }
            let a = d.mul_add(self.weight(flux / d), (-flux).max(0.0));
            return FaceTerm {
                diag: a + flux,
                off: Some((n, a)),
                rhs: 0.0,
                advective: Advective::None,
            };
        }
        let solid = !self.domain.is_fluid(c);
        let half = area * kc / half_c; // A k / (w/2)
        match self.faces[f] {
            ThermalFace::Adiabatic => FaceTerm {
                diag: 0.0,
                off: None,
                rhs: 0.0,
                advective: Advective::None,
            },
            ThermalFace::Temperature(t) => {
                let half = if let Some(previous) = self.conductivity_temperature
                    && let Voxel::Solid(m) = self.domain.voxel_at(c)
                {
                    let material = &self.solids[usize::from(m)];
                    area * material.axis_conductivity()[axis]
                        * material.conductivity_secant_factor(previous[c], t)
                        / half_c
                } else {
                    half
                };
                FaceTerm {
                    diag: half,
                    off: None,
                    rhs: half * t,
                    advective: Advective::None,
                }
            }
            ThermalFace::HeatFlux(q) => FaceTerm {
                diag: 0.0,
                off: None,
                rhs: q * area,
                advective: Advective::None,
            },
            ThermalFace::Convective { h, ambient } => {
                let kc = if let Some(previous) = self.conductivity_temperature
                    && let Voxel::Solid(m) = self.domain.voxel_at(c)
                {
                    let material = &self.solids[usize::from(m)];
                    let axis_k = material.axis_conductivity()[axis];
                    axis_k
                        * material.film_conductivity_factor(previous[c], ambient, axis_k, half_c, h)
                } else {
                    kc
                };
                let u = area / (1.0 / h + half_c / kc);
                FaceTerm {
                    diag: u,
                    off: None,
                    rhs: u * ambient,
                    advective: Advective::None,
                }
            }
            ThermalFace::Inflow { .. } | ThermalFace::Outflow { .. } if solid => FaceTerm {
                diag: 0.0,
                off: None,
                rhs: 0.0,
                advective: Advective::None,
            },
            ThermalFace::Inflow { temperature } => {
                if self.enthalpy.is_some() {
                    let mut term = self.enthalpy_boundary(c, f);
                    term.diag += half;
                    term.rhs += half * temperature;
                    return term;
                }
                let a = half.mul_add(self.weight(flux / half), (-flux).max(0.0));
                let advective = if flux >= 0.0 {
                    Advective::OwnCell(flux)
                } else {
                    Advective::Fixed(flux * temperature)
                };
                FaceTerm {
                    diag: a + flux,
                    off: None,
                    rhs: a * temperature,
                    advective,
                }
            }
            ThermalFace::Outflow {
                backflow_temperature,
            } => {
                if self.enthalpy.is_some() {
                    return self.enthalpy_boundary(c, f);
                }
                if flux >= 0.0 {
                    FaceTerm {
                        diag: flux,
                        off: None,
                        rhs: 0.0,
                        advective: Advective::OwnCell(flux),
                    }
                } else {
                    let j = flux * backflow_temperature;
                    FaceTerm {
                        diag: 0.0,
                        off: None,
                        rhs: -j,
                        advective: Advective::Fixed(j),
                    }
                }
            }
        }
    }
}

/// Solve the steady conjugate energy equation on `domain` with the given
/// fluid, solid material table, flux field, and thermal setup.
///
/// # Errors
/// [`ChtError::InvalidInput`] / [`ChtError::UnknownMaterial`] /
/// [`ChtError::InvalidDomain`] for inadmissible inputs,
/// [`ChtError::FlowThroughClosedFace`] when flow crosses a non-open face,
/// solver refusals, or [`ChtError::Cancelled`].
pub fn solve_energy(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    config: &EnergyConfig,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    solve_energy_inner(domain, fluid, solids, flow, setup, config, None, gate)
}

/// Implicit pseudo-time storage `coefficient_c (T_c - previous_c)` (W) added
/// to every free cell's balance, and the warm start for the Krylov solve.
/// A steady fixed point of the stepped problem is the steady problem; the
/// returned balance excludes the storage term, so it is exact only at that
/// fixed point.
pub(crate) struct PseudoStep<'a> {
    pub coefficient_w_k: &'a [f64],
    pub previous: &'a [f64],
}

/// The energy solve, iterated to a consistent temperature field when a solid
/// declares a temperature-dependent conductivity. Each Picard pass evaluates
/// the integrated bulk laws and their material/contact interface conditions,
/// until the largest temperature change is below 1e-10 of its span.
#[allow(clippy::too_many_arguments)] // physics inputs + step state
pub(crate) fn solve_energy_inner(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    if solids.iter().all(|s| s.conductivity_table.is_empty()) {
        return solve_energy_once(domain, fluid, solids, flow, setup, config, step, None, gate);
    }
    for solid in solids {
        for (i, &(t, k)) in solid.conductivity_table.iter().enumerate() {
            finite("solid.conductivity_table.temperature", t)?;
            finite_positive("solid.conductivity_table.conductivity", k)?;
            if i > 0 && t <= solid.conductivity_table[i - 1].0 {
                return Err(ChtError::InvalidInput {
                    field: "solid.conductivity_table",
                    reason: format!("temperatures must increase ({})", solid.label),
                });
            }
        }
    }
    domain.check_materials(solids.len())?;
    let scales = |temperature: &[f64]| -> Vec<f64> {
        (0..domain.cell_count())
            .map(|c| match domain.voxel_at(c) {
                Voxel::Fluid => 1.0,
                Voxel::Solid(m) => solids[usize::from(m)].conductivity_factor(temperature[c]),
            })
            .collect()
    };
    // First pass: every solid at its table's value for the first declared
    // face, fixed or step temperature (or its reference conductivity).
    let start = match step {
        Some(step) => scales(step.previous),
        None => {
            let reference = setup.faces.iter().find_map(|rule| match *rule {
                ThermalFace::Temperature(t) | ThermalFace::Inflow { temperature: t } => Some(t),
                ThermalFace::Convective { ambient, .. } => Some(ambient),
                _ => None,
            });
            match reference {
                Some(t) => scales(&vec![t; domain.cell_count()]),
                None => vec![1.0; domain.cell_count()],
            }
        }
    };
    let mut solution = solve_energy_once(
        domain,
        fluid,
        solids,
        flow,
        setup,
        config,
        step,
        Some((&start, None)),
        gate,
    )?;
    for _ in 0..200 {
        super::poll(gate)?;
        let next = solve_energy_once(
            domain,
            fluid,
            solids,
            flow,
            setup,
            config,
            step,
            Some((&scales(&solution.temperature), Some(&solution.temperature))),
            gate,
        )?;
        let (lo, hi) = next
            .temperature
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), t| {
                (lo.min(*t), hi.max(*t))
            });
        let change = next
            .temperature
            .iter()
            .zip(&solution.temperature)
            .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
        solution = next;
        if change <= 1e-10 * (hi - lo).max(1e-12) {
            return Ok(solution);
        }
    }
    Err(ChtError::SolverNotConverged {
        system: "energy conductivity",
        iterations: 200,
        relative_residual: f64::NAN,
        tolerance: 1e-10,
    })
}

/// One linear energy solve; `nonlinear` carries per-cell solid conductivity
/// factors and an optional warm start (the previous Picard pass).
#[allow(clippy::too_many_arguments)]
fn solve_energy_once(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>,
    nonlinear: Option<(&[f64], Option<&[f64]>)>,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    solve_energy_field(domain, fluid, solids, flow, setup, config, step, nonlinear, None, None, gate)
}

/// Shared assembly/solve/accounting kernel. An explicit field replaces only
/// molecular conductivity, never material indices, contacts, or heat sources.
#[allow(clippy::too_many_lines, clippy::too_many_arguments)]
fn solve_energy_field(
    domain: &VoxelDomain,
    fluid: &FluidProperties,
    solids: &[SolidMaterial],
    flow: &FlowField,
    setup: &ThermalSetup,
    config: &EnergyConfig,
    step: Option<&PseudoStep<'_>>,
    nonlinear: Option<(&[f64], Option<&[f64]>)>,
    molecular_conductivity: Option<&[[f64; 3]]>,
    enthalpy: Option<AdvectiveEnthalpyField<'_>>,
    gate: &CancelGate,
) -> Result<EnergySolution, ChtError> {
    fluid.validate()?;
    if let Some(field) = enthalpy { field.validate(domain, flow, setup, config, gate)?; }
    for solid in solids {
        finite_positive("solid.conductivity_w_m_k", solid.conductivity_w_m_k)?;
        for k in solid.axis_conductivity() {
            finite_positive("solid.orthotropic_w_m_k", k)?;
        }
    }
    domain.check_materials(solids.len())?;
    if flow.dims() != domain.dims() {
        return Err(ChtError::InvalidDomain {
            reason: format!(
                "flow field {:?} does not match domain {:?}",
                flow.dims(),
                domain.dims()
            ),
        });
    }
    finite_positive("energy.tolerance", config.tolerance)?;
    let cells = domain.cell_count();
    if !(setup.power_w.is_empty() || setup.power_w.len() == cells) {
        return Err(ChtError::InvalidInput {
            field: "thermal.power_w",
            reason: format!("expected 0 or {cells} entries, got {}", setup.power_w.len()),
        });
    }
    for &p in &setup.power_w {
        finite("thermal.power_w", p)?;
    }
    let mut fixed = vec![None; cells];
    for &(c, t) in &setup.fixed_temperature {
        finite("thermal.fixed_temperature", t)?;
        if c >= cells {
            return Err(ChtError::InvalidInput {
                field: "thermal.fixed_temperature",
                reason: format!("cell {c} outside the {cells}-cell domain"),
            });
        }
        fixed[c] = Some(t);
    }
    for (face, rule) in Face3::ALL.into_iter().zip(setup.faces) {
        match rule {
            ThermalFace::Adiabatic => {}
            ThermalFace::Temperature(t) => finite("thermal.face.temperature", t)?,
            ThermalFace::HeatFlux(q) => finite("thermal.face.heat_flux", q)?,
            ThermalFace::Convective { h, ambient } => {
                finite_positive("thermal.face.h", h)?;
                finite("thermal.face.ambient", ambient)?;
            }
            ThermalFace::Inflow { temperature } => finite("thermal.face.inflow", temperature)?,
            ThermalFace::Outflow {
                backflow_temperature,
            } => {
                finite("thermal.face.backflow", backflow_temperature)?;
            }
        }
        if !rule.is_open() {
            let f = face as usize;
            let mut net = 0.0;
            let mut crossed = false;
            for c in 0..cells {
                if domain.neighbor(c, f).is_none() {
                    let q = flow.outward(domain, c, f);
                    crossed |= q != 0.0;
                    net += q;
                }
            }
            if crossed {
                return Err(ChtError::FlowThroughClosedFace {
                    face,
                    net_flux_m3_s: net,
                });
            }
        }
    }
    if !setup.eddy_conductivity_w_m_k.is_empty() {
        if setup.eddy_conductivity_w_m_k.len() != cells {
            return Err(ChtError::InvalidInput {
                field: "thermal.eddy_conductivity_w_m_k",
                reason: format!(
                    "expected {cells} entries, got {}",
                    setup.eddy_conductivity_w_m_k.len()
                ),
            });
        }
        for &k in &setup.eddy_conductivity_w_m_k {
            finite("thermal.eddy_conductivity_w_m_k", k)?;
            if k < 0.0 {
                return Err(ChtError::InvalidInput {
                    field: "thermal.eddy_conductivity_w_m_k",
                    reason: "must be non-negative".into(),
                });
            }
        }
    }
    let eddy = |c: usize| setup.eddy_conductivity_w_m_k.get(c).copied().unwrap_or(0.0);
    let k: Vec<[f64; 3]> = if let Some(field) = molecular_conductivity {
        // Admission is performed by the public field entry point. Keep the
        // existing fluid eddy contribution separate and add it exactly once.
        let mut resolved = field.to_vec();
        for (c, axes) in resolved.iter_mut().enumerate() {
            if domain.is_fluid(c) {
                for axis in axes {
                    *axis += eddy(c);
                    finite_positive("thermal.total_conductivity", *axis)?;
                }
            }
        }
        resolved
    } else {
        (0..cells)
            .map(|c| match domain.voxel_at(c) {
                Voxel::Fluid => [fluid.conductivity_w_m_k + eddy(c); 3],
                Voxel::Solid(m) => {
                    let scale = nonlinear.map_or(1.0, |(scales, _)| scales[c]);
                    solids[usize::from(m)]
                        .axis_conductivity()
                        .map(|k| k * scale)
                }
            })
            .collect()
    };
    let mut contacts = Vec::with_capacity(setup.contacts.len());
    for contact in &setup.contacts {
        let (a, b) = contact.materials;
        if a == b || usize::from(a.max(b)) >= solids.len() {
            return Err(ChtError::InvalidInput {
                field: "thermal.contacts",
                reason: format!(
                    "contact between materials {a} and {b} needs two distinct declared solids"
                ),
            });
        }
        finite(
            "thermal.contacts.resistance_m2_k_w",
            contact.resistance_m2_k_w,
        )?;
        if contact.resistance_m2_k_w < 0.0 {
            return Err(ChtError::InvalidInput {
                field: "thermal.contacts.resistance_m2_k_w",
                reason: "must be non-negative".into(),
            });
        }
        contacts.push(((a.min(b), a.max(b)), contact.resistance_m2_k_w));
    }
    let compact = Compact::new(domain, &k, setup)?;
    let face_contacts = face_contact::prepare(domain, flow, setup, &compact, gate)?;
    let ctx = Context {
        domain,
        flow,
        solids,
        conductivity_temperature: nonlinear.and_then(|(_, temperature)| temperature),
        faces: setup.faces,
        k,
        contacts,
        face_contacts,
        rho_c: fluid.volumetric_heat_capacity(),
        density: fluid.density_kg_m3,
        enthalpy,
        scheme: config.scheme,
    };
    // Validate the actual native or frozen face law, not a different
    // centre-conductivity approximation during nonlinear Picard passes.
    for (index, &(a, b)) in ctx.face_contacts.keys().enumerate() {
        if index.is_multiple_of(256) { super::poll(gate)?; }
        let face = (0..6).find(|&f| domain.neighbor(a, f) == Some(b))
            .expect("admitted contact face");
        if !ctx.diffusive_conductance(a, b, face / 2).is_finite() {
            return Err(ChtError::InvalidInput {
                field: "thermal.face_contacts",
                reason: format!("contact {index}: total face conductance is unrepresentable"),
            });
        }
    }
    let source = |c: usize| {
        if setup.power_w.is_empty() {
            0.0
        } else {
            setup.power_w[c]
        }
    };
    let parts = setup.compact_components.len();
    let unknowns = cells + parts;

    let mut coo = Coo::new(unknowns, unknowns);
    let mut b = vec![0.0f64; unknowns];
    let mut max_cell_peclet = 0.0f64;
    // The box grid is face-connected, so one anchor anywhere makes the
    // operator nonsingular; without one the problem is pure Neumann.
    let mut anchored = !setup.fixed_temperature.is_empty();
    // Sinks per cell, validated once.
    let mut sinks: Vec<(f64, f64)> = vec![(0.0, 0.0); cells];
    for sink in &setup.cell_sinks {
        if sink.cell >= cells {
            return Err(ChtError::InvalidInput {
                field: "thermal.cell_sinks",
                reason: format!("cell {} is outside the domain", sink.cell),
            });
        }
        finite("thermal.cell_sinks.conductance_w_k", sink.conductance_w_k)?;
        finite("thermal.cell_sinks.temperature_k", sink.temperature_k)?;
        if sink.conductance_w_k < 0.0 {
            return Err(ChtError::InvalidInput {
                field: "thermal.cell_sinks.conductance_w_k",
                reason: "must be non-negative".into(),
            });
        }
        let slot = &mut sinks[sink.cell];
        slot.0 += sink.conductance_w_k;
        slot.1 = sink.conductance_w_k.mul_add(sink.temperature_k, slot.1);
        anchored |= sink.conductance_w_k > 0.0;
    }
    for c in 0..cells {
        if c.is_multiple_of(4096) {
            super::poll(gate)?;
        }
        if let Some(j) = compact.owner(c) {
            // Collapsed component cell: it reports its junction temperature.
            coo.push(c, c, 1.0);
            coo.push(c, cells + j, -1.0);
            continue;
        }
        if let Some(t) = fixed[c] {
            coo.push(c, c, 1.0);
            b[c] = t;
            continue;
        }
        let (mut diag, mut rhs) = match step {
            Some(step) => {
                let coefficient = step.coefficient_w_k[c];
                finite_positive("thermal.storage.coefficient_w_k", coefficient)?;
                finite("thermal.storage.previous_k", step.previous[c])?;
                // Physical storage fixes the absolute temperature even when
                // every boundary is insulated or supplies only a heat flux.
                anchored = true;
                (coefficient, coefficient.mul_add(step.previous[c], source(c)))
            }
            None => (0.0, source(c)),
        };
        // Cell sinks: G (T - T_sink) leaves the cell.
        diag += sinks[c].0;
        rhs += sinks[c].1;
        // Junction resistors reaching this cell.
        for &(j, g) in compact.links_of(c) {
            diag += g;
            coo.push(c, cells + j, -g);
        }
        for f in 0..6 {
            if compact.toward_component(c, f) {
                continue;
            }
            let term = ctx.checked_term(c, f)?;
            diag += term.diag;
            rhs += term.rhs;
            anchored |= term.off.is_none() && term.diag > 0.0;
            if let Some((n, a)) = term.off {
                coo.push(c, n, -a);
                if domain.is_fluid(c) && domain.is_fluid(n) {
                    // A k / (centre distance).
                    let axis = f / 2;
                    let d = domain.face_area(c, axis) * ctx.k[c][axis]
                        / (0.5 * (domain.widths(c)[axis] + domain.widths(n)[axis]));
                    let q = flow.outward(domain, c, f);
                    let rho_c = ctx.enthalpy.map_or(ctx.rho_c, |field| {
                        ctx.density * field.cells[if q >= 0.0 { c } else { n }]
                            .expect("admitted fluid donor").heat_capacity_j_kg_k
                    });
                    let flux = rho_c * q;
                    max_cell_peclet = max_cell_peclet.max(flux.abs() / d);
                }
            }
        }
        if diag <= 0.0 {
            return Err(ChtError::InvalidInput {
                field: "thermal.faces",
                reason: format!(
                    "cell {c} has no temperature anchor (isolated adiabatic region or inflow-only stencil)"
                ),
            });
        }
        coo.push(c, c, diag);
        b[c] = rhs;
    }
    for (j, part) in setup.compact_components.iter().enumerate() {
        let row = cells + j;
        let mut diag = 0.0;
        for &(n, g, _) in &compact.links[j] {
            diag += g;
            coo.push(row, n, -g);
        }
        coo.push(row, row, diag);
        b[row] = part.power_w;
    }
    if !anchored {
        return Err(ChtError::InvalidInput {
            field: "thermal.faces",
            reason: "no temperature anchor: declare a Temperature, Convective, Inflow or flowing Outflow face, or a fixed cell".into(),
        });
    }
    let a = scale_rows(&coo, &mut b);
    let nonzeros = a.nnz();
    let guess = setup
        .faces
        .iter()
        .find_map(|rule| match *rule {
            ThermalFace::Inflow { temperature } => Some(temperature),
            ThermalFace::Temperature(t) => Some(t),
            ThermalFace::Convective { ambient, .. } => Some(ambient),
            _ => None,
        })
        .or_else(|| setup.fixed_temperature.first().map(|&(_, t)| t))
        .unwrap_or(0.0);
    let warm = nonlinear.and_then(|(_, warm)| warm);
    let mut temperature = match (warm, step) {
        (Some(previous), _) => {
            let mut start = previous.to_vec();
            start.extend(
                setup
                    .compact_components
                    .iter()
                    .map(|part| previous[domain.index(part.lo[0], part.lo[1], part.lo[2])]),
            );
            start
        }
        (None, Some(step)) => {
            let mut warm = step.previous.to_vec();
            warm.extend(
                setup
                    .compact_components
                    .iter()
                    .map(|part| step.previous[domain.index(part.lo[0], part.lo[1], part.lo[2])]),
            );
            warm
        }
        (None, None) => vec![guess; unknowns],
    };
    let outcome = bicgstab_ilu0(
        "energy",
        &a,
        &b,
        &mut temperature,
        config.tolerance,
        config.max_iterations,
        gate,
    )?;

    // Independent accounting from boundary and fixed-cell fluxes.
    let (mut source_w, mut fixed_w, mut out_w, mut adv_w, mut scale) = (0.0, 0.0, 0.0, 0.0, 0.0f64);
    let flux_of = |term: &FaceTerm, c: usize| {
        let off = term.off.map_or(0.0, |(n, a)| a * temperature[n]);
        term.diag.mul_add(temperature[c], -off) - term.rhs
    };
    let mut sink_w = 0.0f64;
    let mut solid_fluid_links = Vec::new();
    for c in 0..cells {
        if c.is_multiple_of(4096) { super::poll(gate)?; }
        if !domain.is_fluid(c) && compact.owner(c).is_none() {
            for f in 0..6 {
                if let Some(n) = domain.neighbor(c, f) && domain.is_fluid(n) {
                    let term = ctx.checked_term(c, f)?;
                    let (_, g) = term.off.expect("interior solid/fluid face");
                    solid_fluid_links.push((c, n, g));
                }
            }
        }
        let sink = sinks[c].0.mul_add(temperature[c], -sinks[c].1);
        sink_w += sink;
        scale += sink.abs();
        if fixed[c].is_some() {
            // A fixed cell also feeds its own sinks and junction links.
            fixed_w += sink;
            for f in 0..6 {
                if compact.toward_component(c, f) {
                    continue;
                }
                let term = ctx.checked_term(c, f)?;
                fixed_w += flux_of(&term, c);
            }
            for &(j, g) in compact.links_of(c) {
                fixed_w += g * (temperature[c] - temperature[cells + j]);
            }
        } else {
            source_w += source(c);
            scale += source(c).abs();
        }
        for f in 0..6 {
            if domain.neighbor(c, f).is_some() {
                continue;
            }
            let term = ctx.checked_term(c, f)?;
            let advective = match term.advective {
                Advective::None => 0.0,
                Advective::OwnCell(q) => q * temperature[c],
                Advective::Affine { mass_flow_kg_s, enthalpy } =>
                    mass_flow_kg_s * enthalpy.at_temperature(temperature[c]),
                Advective::Fixed(j) => j,
            };
            // Centered enthalpy/diffusion avoids cancelling large affine
            // matrix coefficients when auditing a small physical heat flux.
            let j = if ctx.enthalpy.is_some() && domain.is_fluid(c) {
                match setup.faces[f] {
                    ThermalFace::Inflow { temperature: inlet } => advective
                        + domain.face_area(c, f / 2) * ctx.k[c][f / 2]
                            / (0.5 * domain.widths(c)[f / 2]) * (temperature[c] - inlet),
                    ThermalFace::Outflow { .. } => advective,
                    _ => flux_of(&term, c),
                }
            } else { flux_of(&term, c) };
            if fixed[c].is_none() {
                out_w += j;
            } else {
                // The fixed cell's boundary flux is part of its own
                // injection; only its interior share entered the domain.
                fixed_w -= j;
            }
            scale += j.abs();
            adv_w += advective;
        }
    }
    for part in &setup.compact_components {
        source_w += part.power_w;
        scale += part.power_w.abs();
    }
    scale += fixed_w.abs();
    let residual_w = source_w + fixed_w - out_w - sink_w;
    let balance = EnergyBalance {
        source_w,
        fixed_cell_injection_w: fixed_w,
        boundary_outflow_w: out_w,
        advective_outflow_w: adv_w,
        sink_outflow_w: sink_w,
        residual_w,
        relative_residual: if scale > 0.0 {
            residual_w.abs() / scale
        } else {
            0.0
        },
    };
    let junctions = setup
        .compact_components
        .iter()
        .enumerate()
        .map(|(j, part)| {
            let t = temperature[cells + j];
            let (mut case_w, mut board_w) = (0.0, 0.0);
            for &(n, g, board) in &compact.links[j] {
                let q = g * (t - temperature[n]);
                if board {
                    board_w += q;
                } else {
                    case_w += q;
                }
            }
            JunctionSolution {
                temperature_k: t,
                case_w,
                board_w,
                cells: (part.lo, part.hi),
                links: compact.links[j].iter().map(|&(n, g, _)| (n, g)).collect(),
            }
        })
        .collect();
    temperature.truncate(cells);
    let report = EnergyReport {
        unknowns,
        nonzeros,
        iterations: outcome.iterations,
        relative_residual: outcome.relative_residual,
        max_cell_peclet,
        scheme: config.scheme,
        balance,
    };
    Ok(EnergySolution {
        temperature,
        conductivity: ctx.k,
        face_contact_resistance_m2_k_w: ctx.face_contacts,
        solid_fluid_links,
        junctions,
        report,
    })
}

/// Compact components prepared for assembly: the owning component of each
/// cell, each junction's links `(outside cell, conductance, board side)`,
/// and the same links indexed by outside cell.
struct Compact {
    owner: Vec<usize>,
    links: Vec<Vec<(usize, f64, bool)>>,
    by_cell: std::collections::BTreeMap<usize, Vec<(usize, f64)>>,
    dims: [usize; 3],
}

impl Compact {
    fn new(domain: &VoxelDomain, k: &[[f64; 3]], setup: &ThermalSetup) -> Result<Self, ChtError> {
        let parts = &setup.compact_components;
        let cells = domain.cell_count();
        let n = domain.dims();
        let refuse = |reason: String| ChtError::InvalidInput {
            field: "thermal.compact_components",
            reason,
        };
        let mut owner = if parts.is_empty() {
            Vec::new()
        } else {
            vec![usize::MAX; cells]
        };
        for (i, part) in parts.iter().enumerate() {
            if (0..3).any(|a| part.lo[a] >= part.hi[a] || part.hi[a] > n[a]) {
                return Err(refuse(format!(
                    "component {i}: cells {:?}..{:?} must be non-empty within {n:?}",
                    part.lo, part.hi
                )));
            }
            finite("thermal.compact_components.power_w", part.power_w)?;
            finite_positive(
                "thermal.compact_components.junction_to_case_k_w",
                part.junction_to_case_k_w,
            )?;
            finite_positive(
                "thermal.compact_components.junction_to_board_k_w",
                part.junction_to_board_k_w,
            )?;
            for z in part.lo[2]..part.hi[2] {
                for y in part.lo[1]..part.hi[1] {
                    for x in part.lo[0]..part.hi[0] {
                        let c = domain.index(x, y, z);
                        if domain.is_fluid(c) {
                            return Err(refuse(format!(
                                "component {i}: cell ({x}, {y}, {z}) is fluid; the box must be solid"
                            )));
                        }
                        if owner[c] != usize::MAX {
                            return Err(refuse(format!("components {} and {i} overlap", owner[c])));
                        }
                        owner[c] = i;
                    }
                }
            }
        }
        let mut links = Vec::with_capacity(parts.len());
        let mut by_cell: std::collections::BTreeMap<usize, Vec<(usize, f64)>> =
            std::collections::BTreeMap::new();
        for (i, part) in parts.iter().enumerate() {
            let board = part.board_face as usize;
            let axis = board / 2;
            let area = part.face_area(domain);
            let mut list = Vec::new();
            for z in part.lo[2]..part.hi[2] {
                for y in part.lo[1]..part.hi[1] {
                    for x in part.lo[0]..part.hi[0] {
                        let c = domain.index(x, y, z);
                        for f in [2 * axis, 2 * axis + 1] {
                            let resistance = if f == board {
                                part.junction_to_board_k_w
                            } else {
                                part.junction_to_case_k_w
                            };
                            match domain.neighbor(c, f) {
                                Some(m) if owner[m] == i => {}
                                Some(m) if owner[m] != usize::MAX => {
                                    return Err(refuse(format!(
                                        "component {i} sits directly on component {}",
                                        owner[m]
                                    )));
                                }
                                None => {
                                    return Err(refuse(format!(
                                        "component {i}: its case or board face lies on the domain boundary"
                                    )));
                                }
                                Some(m) => {
                                    // The resistor's area share in series
                                    // with the outside cell's half width.
                                    let face_area = domain.face_area(m, axis);
                                    let g = 1.0
                                        / (resistance * area / face_area
                                            + 0.5 * domain.widths(m)[axis]
                                                / (k[m][axis] * face_area));
                                    list.push((m, g, f == board));
                                    by_cell.entry(m).or_default().push((i, g));
                                }
                            }
                        }
                    }
                }
            }
            links.push(list);
        }
        if !owner.is_empty() {
            for &(c, _) in &setup.fixed_temperature {
                if owner.get(c).is_some_and(|&o| o != usize::MAX) {
                    return Err(refuse(format!("cell {c} is fixed inside a component")));
                }
            }
            for sink in &setup.cell_sinks {
                if owner.get(sink.cell).is_some_and(|&o| o != usize::MAX) {
                    return Err(refuse(format!(
                        "cell {} carries a sink inside a component",
                        sink.cell
                    )));
                }
            }
            for (c, &p) in setup.power_w.iter().enumerate() {
                if p != 0.0 && owner[c] != usize::MAX {
                    return Err(refuse(format!(
                        "cell {c} carries power inside a component; declare it as the component's power_w"
                    )));
                }
            }
        }
        Ok(Self {
            owner,
            links,
            by_cell,
            dims: n,
        })
    }

    fn owner(&self, c: usize) -> Option<usize> {
        self.owner.get(c).copied().filter(|&o| o != usize::MAX)
    }

    /// Whether face `f` of cell `c` leads into a component cell (those faces
    /// carry no ordinary conduction; the junction links replace them).
    fn toward_component(&self, c: usize, f: usize) -> bool {
        if self.owner.is_empty() {
            return false;
        }
        let mut at = [
            c % self.dims[0],
            (c / self.dims[0]) % self.dims[1],
            c / (self.dims[0] * self.dims[1]),
        ];
        let axis = f / 2;
        if f % 2 == 1 {
            if at[axis] + 1 >= self.dims[axis] {
                return false;
            }
            at[axis] += 1;
        } else {
            if at[axis] == 0 {
                return false;
            }
            at[axis] -= 1;
        }
        let n = (at[2] * self.dims[1] + at[1]) * self.dims[0] + at[0];
        self.owner[n] != usize::MAX
    }

    fn links_of(&self, c: usize) -> &[(usize, f64)] {
        self.by_cell.get(&c).map_or(&[], Vec::as_slice)
    }
}
