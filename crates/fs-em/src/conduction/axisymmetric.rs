//! Axisymmetric Ohmic transport on the original meridional `(r,z)` complex.
//!
//! The native conductivity Hodge integrates `dV = angular_span * r dr dz`.
//! Voltage, inward current traces and cellwise `E.J` retain the same source
//! admission, checked scalar solve and independent current/work gates as the
//! planar owner. Surface loads integrate the radial P1 moments exactly.
//!
//! Fields and material axes are local `[radial, axial]` components, independent
//! of azimuth. An oriented meridional law must have a decoupled azimuthal
//! direction; this is not a revolution of fixed Cartesian anisotropy. Current,
//! heat and resistance refer to the supplied angular span (TAU for a full body).
//! A smaller span is a symmetry sector, not independently driven wedge faces.
//! The zero-area axis has natural regularity and cannot be an electrode.
//!
//! This is stationary fixed-geometry conduction. It does not add azimuthal
//! current, induction, Hall response, moving contacts or thermal evolution.
//! Material changes require a newly prepared source-bound operator. The private
//! scalar core cannot escape through the planar motion or contact APIs.

use super::{
    ConductionControls, ConductionCurrentTrace, ConductionError, ConductionField,
    ConductionPowerSource, PlanarConduction, PlanarResistance, poll, scalar_potential,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::conductor::{
    ResolvedConductorState,
    planar::{PlanarConductivity, PlanarOhmicResponse},
};
use fs_phs::circuit::Branch;
use fs_rep_mesh::TriComplex2;
use fs_solver::SolveReport;

mod matrix;
pub use matrix::{
    AxisymmetricConductanceControls, AxisymmetricConductanceMatrix, AxisymmetricMixedConductance,
    AxisymmetricMixedConductanceControls,
};

type Result<T> = core::result::Result<T, ConductionError>;

/// Original radial geometry and complete, immutable conductivity sources.
pub struct AxisymmetricConduction<C = ResolvedConductorState> {
    core: PlanarConduction<C>,
}

impl<C: PlanarOhmicResponse> AxisymmetricConduction<C> {
    /// Prepare scalar or explicitly oriented meridional conductivity in S/m.
    /// Negative radii and invalid physical measures are refused by TriComplex2.
    /// All oriented sources must name one common meridional geometry frame.
    pub fn new(
        mesh: TriComplex2,
        conductors: Vec<C>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        if !mesh.metric().is_axisymmetric() {
            return Err(ConductionError::Invalid(
                "axisymmetric conduction requires an explicit radial metric",
            ));
        }
        Ok(Self {
            core: PlanarConduction::assemble_metric(
                mesh,
                conductors,
                limits,
                cancelled,
                |source| source.has_admitted_ohmic_state(),
            )?,
        })
    }

    /// Metre-valued `[radius, axial_coordinate]` vertices and swept measures.
    pub fn mesh(&self) -> &TriComplex2 {
        self.core.mesh()
    }

    /// Complete retained source/state snapshots, in original triangle order.
    pub fn conductors(&self) -> &[C] {
        self.core.conductors()
    }

    /// Scalar or full meridional conductivity, without a tensor surrogate.
    pub fn conductivity_coefficients(&self) -> &[PlanarConductivity] {
        self.core.conductivity_coefficients()
    }

    /// Exterior voltage traces; all unlisted physical surfaces are insulating.
    /// Every connected component needs a positive-area prescribed trace.
    pub fn solve(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'_, C>> {
        self.solve_with_current_traces(boundary_voltage_v, &[], controls, cancelled)
    }

    /// Uniform inward current density [A/m2] on each named surface of
    /// revolution. The native field determines free voltages. These are
    /// density traces, not equipotential total-current electrodes.
    pub fn solve_with_current_traces(
        &self,
        boundary_voltage_v: &[(usize, f64)],
        current_traces: &[ConductionCurrentTrace],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'_, C>> {
        let field = self.core.solve_with_current_traces(
            boundary_voltage_v,
            current_traces,
            controls,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self,
            field,
        })
    }

    /// Native resistance between two disjoint complete physical electrodes.
    /// Use the returned response for voltage, total-current or power drives.
    pub fn resistance(
        &self,
        positive_terminal: &[usize],
        negative_terminal: &[usize],
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricResistance<'_, C>> {
        self.resistance_at_reference_voltage(
            positive_terminal,
            negative_terminal,
            1.0,
            controls,
            cancelled,
        )
    }

    /// The same fixed source/geometry response at a positive probe voltage
    /// inside the material domain. Preparation contributes no separate heat.
    pub fn resistance_at_reference_voltage(
        &self,
        positive_terminal: &[usize],
        negative_terminal: &[usize],
        reference_voltage_v: f64,
        controls: ConductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricResistance<'_, C>> {
        admit_terminal_traces(
            self.mesh(),
            &[positive_terminal, negative_terminal],
            cancelled,
        )?;
        let response = self.core.resistance_at_reference_voltage(
            positive_terminal,
            negative_terminal,
            reference_voltage_v,
            controls,
            cancelled,
        )?;
        Ok(AxisymmetricResistance {
            problem: self,
            response,
        })
    }
}

// Validate every electrode separately: an axis-only terminal must not borrow
// positive area from the neighboring edge of another electrode.
fn admit_terminal_traces(
    mesh: &TriComplex2,
    terminals: &[&[usize]],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    if terminals.len() < 2 || terminals.len() > mesh.vertices().len() {
        return Err(ConductionError::Invalid(
            "bounded physical axisymmetric terminal traces required",
        ));
    }
    // Native admission checks every index, duplicate, overlap and edge first.
    scalar_potential::terminal_boundary(mesh, terminals, &vec![0.0; terminals.len()], cancelled)?;
    let mut owners = vec![None; mesh.vertices().len()];
    for (terminal, vertices) in terminals.iter().enumerate() {
        for &vertex in *vertices {
            poll(cancelled)?;
            owners[vertex] = Some(terminal);
        }
    }
    admit_physical_traces(mesh, &owners, cancelled)
}

// Integral(N_i r ds)/integral(r ds) = (1+r_i/(r_0+r_1))/3.
// Scale radii before summing; the admitted edge has positive physical area.
pub(super) fn surface_current_loads(radii: [f64; 2], current_a: f64) -> [f64; 2] {
    let scale = radii[0].max(radii[1]);
    let scaled = radii.map(|r| r / scale);
    let sum = scaled[0] + scaled[1];
    scaled.map(|r| current_a * ((1.0 + r / sum) / 3.0))
}

pub(super) fn admit_prescribed_traces(
    mesh: &TriComplex2,
    prescribed: &[Option<f64>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let mut owners = Vec::with_capacity(prescribed.len());
    for value in prescribed {
        poll(cancelled)?;
        owners.push(value.map(|_| 0));
    }
    admit_physical_traces(mesh, &owners, cancelled)
}

// An axis line has positive meridional length but zero physical surface area.
// A terminal cannot borrow area from an adjacent DIFFERENT terminal.
fn admit_physical_traces(
    mesh: &TriComplex2,
    owners: &[Option<usize>],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    poll(cancelled)?;
    let trace = mesh
        .boundary_trace()
        .map_err(|_| ConductionError::Invalid("invalid axisymmetric electrical boundary trace"))?;
    let mut supported = vec![false; owners.len()];
    for edge in &trace.edges {
        poll(cancelled)?;
        let [a, b] = edge.oriented_vertices.map(|v| v as usize);
        if owners[a].is_some()
            && owners[a] == owners[b]
            && mesh
                .edge_measure(edge.global_edge)
                .is_some_and(|area| area > 0.0)
        {
            supported[a] = true;
            supported[b] = true;
        }
    }
    if owners
        .iter()
        .zip(supported)
        .any(|(owner, ok)| owner.is_some() && !ok)
    {
        return Err(ConductionError::Invalid(
            "axisymmetric electrodes require physical surface area, not the symmetry axis",
        ));
    }
    poll(cancelled)
}

/// Accepted meridional current and Joule heat tied to original source geometry.
pub struct AxisymmetricConductionField<'a, C = ResolvedConductorState> {
    problem: &'a AxisymmetricConduction<C>,
    field: ConductionField<'a, C>,
}

impl<C: PlanarOhmicResponse> AxisymmetricConductionField<'_, C> {
    pub fn problem(&self) -> &AxisymmetricConduction<C> {
        self.problem
    }

    pub fn voltage_v(&self) -> &[f64] {
        self.field.voltage_v()
    }

    /// Physical recovery uses these component-relative potentials, retaining
    /// drops that may round away in the absolute voltage reporting view.
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.field.relative_voltage_v()
    }

    pub fn reference_voltage_v(&self) -> &[f64] {
        self.field.reference_voltage_v()
    }

    /// Original cell `[E_r, E_z]` in V/m, with no azimuthal component.
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        self.field.electric_field_v_m()
    }

    /// Original cell `[J_r, J_z]` in A/m2, without interface smoothing.
    pub fn current_density_a_m2(&self) -> &[[f64; 2]] {
        self.field.current_density_a_m2()
    }

    /// Integral(E.J dV) on each original swept cell [W]. This is the same
    /// irreversible loss as the native resistor, not a second thermal credit.
    pub fn cell_power_w(&self) -> &[f64] {
        self.field.cell_power_w()
    }

    pub fn power_w(&self) -> f64 {
        self.field.power_w()
    }

    /// Combined voltage reaction and prescribed trace current, positive inward.
    pub fn nodal_boundary_current_a(&self) -> &[Option<f64>] {
        self.field.nodal_boundary_current_a()
    }

    /// Exact radial P1 moments of the prescribed density traces [A].
    pub fn nodal_applied_current_a(&self) -> &[Option<f64>] {
        self.field.nodal_applied_current_a()
    }

    pub fn boundary_current_a(&self) -> f64 {
        self.field.boundary_current_a()
    }

    pub fn maximum_free_current_residual_a(&self) -> f64 {
        self.field.maximum_free_current_residual_a()
    }

    /// Joule power minus supplied voltage/current work [W].
    pub fn power_balance_residual_w(&self) -> f64 {
        self.field.power_balance_residual_w()
    }

    pub fn report(&self) -> Option<&SolveReport> {
        self.field.report()
    }
}

/// One native reference solve; repeated drives recover fields without Krylov.
/// Reprepare after geometry, temperature or another material-state change.
pub struct AxisymmetricResistance<'a, C = ResolvedConductorState> {
    problem: &'a AxisymmetricConduction<C>,
    response: PlanarResistance<'a, C>,
}

impl<'a, C: PlanarOhmicResponse> AxisymmetricResistance<'a, C> {
    pub fn resistance_ohm(&self) -> f64 {
        self.response.resistance_ohm()
    }

    pub fn reference_voltage_v(&self) -> f64 {
        self.response.reference_voltage_v()
    }

    /// The same dissipative component; do not count its field heat twice.
    /// A circuit owner must admit each actual drive with `field_at_voltage`.
    pub fn circuit_branch(&self) -> Branch {
        self.response.circuit_branch()
    }

    pub fn field_at_voltage(
        &self,
        voltage_v: f64,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'a, C>> {
        let field = self.response.field_at_voltage(
            voltage_v,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }

    /// Signed equipotential total current, with both reactions and I*V checked.
    pub fn field_at_current(
        &self,
        current_a: f64,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'a, C>> {
        let field = self.response.field_at_current(
            current_a,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }

    /// Ideal regulated power with hard voltage/current compliance limits.
    pub fn field_at_power(
        &self,
        source: ConductionPowerSource,
        current_tolerance_a: f64,
        power_tolerance_w: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricConductionField<'a, C>> {
        let field = self.response.field_at_power(
            source,
            current_tolerance_a,
            power_tolerance_w,
            cancelled,
        )?;
        Ok(AxisymmetricConductionField {
            problem: self.problem,
            field,
        })
    }
}
