//! Steady conservative finite-volume energy equation over fluid and solid
//! voxels.
//!
//! Every face contributes one outward total energy flux written as
//! `J = d T_P - a T_N - r` (Patankar, *Numerical Heat Transfer and Fluid
//! Flow*, 1980, §5.2–5.4): for an interior face with outward heat-capacity
//! flux `F = rho c_p Q` and diffusive conductance `D = A k_h / dx`,
//! `a = D A(|F|/D) + max(-F, 0)` and `d = a + F`; the two sides of one face
//! therefore carry equal and opposite `J` for ANY temperatures, so the global
//! energy balance telescopes exactly. `k_h = 2 k_P k_N / (k_P + k_N)` is the
//! series (harmonic) conductance of two half cells, exact for
//! piecewise-constant conductivity meeting at the face.

use fs_exec::CancelGate;
use fs_sparse::Coo;

use super::domain::{FluidProperties, SolidMaterial, Voxel, VoxelDomain};
use super::flow::{FlowField, scale_rows};
use super::krylov::bicgstab_ilu0;
use super::{ChtError, finite, finite_positive};
use crate::d3q19::Face3;

/// Thermal rule on one domain face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalFace {
    /// Zero heat flux. Refuses if flow crosses the face.
    Adiabatic,
    /// Prescribed temperature (K) on the face plane, half a cell from the
    /// boundary cell centre.
    Temperature(f64),
    /// Prescribed heat flux INTO the domain, W/m^2.
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
    fn contains(&self, at: [usize; 3]) -> bool {
        (0..3).all(|a| at[a] >= self.lo[a] && at[a] < self.hi[a])
    }

    /// Area of the case-top (= board) face, m^2.
    fn face_area(&self, dx: f64) -> f64 {
        let axis = self.board_face as usize / 2;
        (0..3)
            .filter(|&a| a != axis)
            .map(|a| (self.hi[a] - self.lo[a]) as f64 * dx)
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
            cell_sinks: Vec::new(),
            eddy_conductivity_w_m_k: Vec::new(),
            compact_components: Vec::new(),
        }
    }

    /// Spread `total_w` uniformly over the cells whose centre satisfies
    /// `inside`; returns the number of heated cells.
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
        if !cells.is_empty() {
            let each = total_w / cells.len() as f64;
            for c in &cells {
                self.power_w[*c] += each;
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
    /// enthalpy referenced to 0 K, which is exact because net boundary
    /// mass flow vanishes).
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
    /// Conductivity per cell along x, y, z used by the solve, W/(m K).
    pub conductivity: Vec<[f64; 3]>,
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
    /// fluid/solid face, W.
    #[must_use]
    pub fn solid_to_fluid_heat_w(&self, domain: &VoxelDomain) -> f64 {
        let dx = domain.dx();
        let mut total = 0.0;
        // A collapsed component cell holds its junction's temperature and
        // reaches its neighbours only through the junction's links.
        let collapsed = |c: usize| {
            let at = domain.coords(c);
            self.junctions
                .iter()
                .any(|j| (0..3).all(|a| at[a] >= j.cells.0[a] && at[a] < j.cells.1[a]))
        };
        for junction in &self.junctions {
            for &(n, g) in &junction.links {
                if domain.is_fluid(n) {
                    total += g * (junction.temperature_k - self.temperature[n]);
                }
            }
        }
        for c in 0..domain.cell_count() {
            if domain.is_fluid(c) || collapsed(c) {
                continue;
            }
            for f in 0..6 {
                if let Some(n) = domain.neighbor(c, f)
                    && domain.is_fluid(n)
                {
                    let axis = f / 2;
                    let (kc, kn) = (self.conductivity[c][axis], self.conductivity[n][axis]);
                    let d = dx * 2.0 * kc * kn / (kc + kn);
                    total += d * (self.temperature[c] - self.temperature[n]);
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
    Fixed(f64),
}

struct Context<'a> {
    domain: &'a VoxelDomain,
    flow: &'a FlowField,
    faces: [ThermalFace; 6],
    /// Conductivity along x, y, z per cell.
    k: Vec<[f64; 3]>,
    /// Contact resistances keyed by the ordered material pair.
    contacts: Vec<((u16, u16), f64)>,
    rho_c: f64,
    scheme: ConvectionScheme,
}

impl Context<'_> {
    fn contact(&self, c: usize, n: usize) -> f64 {
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

    fn term(&self, c: usize, f: usize) -> FaceTerm {
        let dx = self.domain.dx();
        let area = dx * dx;
        let flux = self.rho_c * self.flow.outward(self.domain, c, f);
        let kc = self.k[c][f / 2];
        if let Some(n) = self.domain.neighbor(c, f) {
            let kn = self.k[n][f / 2];
            // Series resistance of the two half cells and any contact
            // between their materials: exact for piecewise-constant k.
            let contact = self.contact(c, n);
            let d = area / (0.5 * dx / kc + contact + 0.5 * dx / kn);
            let a = d.mul_add(self.weight(flux / d), (-flux).max(0.0));
            return FaceTerm {
                diag: a + flux,
                off: Some((n, a)),
                rhs: 0.0,
                advective: Advective::None,
            };
        }
        let solid = !self.domain.is_fluid(c);
        let half = 2.0 * dx * kc; // A k / (dx/2)
        match self.faces[f] {
            ThermalFace::Adiabatic => FaceTerm {
                diag: 0.0,
                off: None,
                rhs: 0.0,
                advective: Advective::None,
            },
            ThermalFace::Temperature(t) => FaceTerm {
                diag: half,
                off: None,
                rhs: half * t,
                advective: Advective::None,
            },
            ThermalFace::HeatFlux(q) => FaceTerm {
                diag: 0.0,
                off: None,
                rhs: q * area,
                advective: Advective::None,
            },
            ThermalFace::Convective { h, ambient } => {
                let u = area / (1.0 / h + 0.5 * dx / kc);
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

#[allow(clippy::too_many_lines, clippy::too_many_arguments)] // admission, assembly, solve, balance
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
    fluid.validate()?;
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
    let k: Vec<[f64; 3]> = (0..cells)
        .map(|c| match domain.voxel_at(c) {
            Voxel::Fluid => [fluid.conductivity_w_m_k + eddy(c); 3],
            Voxel::Solid(m) => solids[usize::from(m)].axis_conductivity(),
        })
        .collect();
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
    let ctx = Context {
        domain,
        flow,
        faces: setup.faces,
        k,
        contacts,
        rho_c: fluid.volumetric_heat_capacity(),
        scheme: config.scheme,
    };
    let source = |c: usize| {
        if setup.power_w.is_empty() {
            0.0
        } else {
            setup.power_w[c]
        }
    };
    let compact = Compact::new(domain, &ctx.k, setup)?;
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
            Some(step) => (
                step.coefficient_w_k[c],
                step.coefficient_w_k[c].mul_add(step.previous[c], source(c)),
            ),
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
            let term = ctx.term(c, f);
            diag += term.diag;
            rhs += term.rhs;
            anchored |= term.off.is_none() && term.diag > 0.0;
            if let Some((n, a)) = term.off {
                coo.push(c, n, -a);
                if domain.is_fluid(c) && domain.is_fluid(n) {
                    let d = domain.dx() * ctx.k[c][f / 2];
                    let flux = ctx.rho_c * flow.outward(domain, c, f);
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
    let mut temperature = match step {
        Some(step) => {
            let mut warm = step.previous.to_vec();
            warm.extend(
                setup
                    .compact_components
                    .iter()
                    .map(|part| step.previous[domain.index(part.lo[0], part.lo[1], part.lo[2])]),
            );
            warm
        }
        None => vec![guess; unknowns],
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
    for c in 0..cells {
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
                let term = ctx.term(c, f);
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
            let term = ctx.term(c, f);
            let j = flux_of(&term, c);
            if fixed[c].is_none() {
                out_w += j;
            } else {
                // The fixed cell's boundary flux is part of its own
                // injection; only its interior share entered the domain.
                fixed_w -= j;
            }
            scale += j.abs();
            adv_w += match term.advective {
                Advective::None => 0.0,
                Advective::OwnCell(q) => q * temperature[c],
                Advective::Fixed(j) => j,
            };
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
        let dx = domain.dx();
        let face_area = dx * dx;
        let mut links = Vec::with_capacity(parts.len());
        let mut by_cell: std::collections::BTreeMap<usize, Vec<(usize, f64)>> =
            std::collections::BTreeMap::new();
        for (i, part) in parts.iter().enumerate() {
            let board = part.board_face as usize;
            let axis = board / 2;
            let area = part.face_area(dx);
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
                                    let g = 1.0
                                        / (resistance * area / face_area
                                            + 0.5 * dx / (k[m][axis] * face_area));
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
