//! Conserved cell mass and mass*h on the original projected volume fluxes.

use super::{
    CancelGate, ChtError, DensityProjection, DensityProjectionControls, FlowFace, FlowField,
    VoxelDomain, finite, finite_positive, poll, poll_every, scale_rows,
};
use crate::conjugate::krylov::bicgstab_ilu0;
use fs_sparse::{Coo, Csr};

mod momentum;
pub use momentum::{
    MaterialMomentumControls, MaterialMomentumInput, MaterialMomentumReport, MaterialMomentumStep,
};

/// Incoming material in the same physical and caloric convention as the cells.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdvectedMaterial {
    /// Positive physical density [kg/m³].
    pub density_kg_m3: f64,
    /// Specific enthalpy [J/kg], including any already absorbed latent heat.
    pub specific_enthalpy_j_kg: f64,
}

/// Explicit material accompanying a pressure projection's volume source.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum MaterialVolumeSource {
    /// No local material source.
    #[default]
    None,
    /// Positive supplied volume rate [m³/s] with an authored incoming state.
    Inflow {
        volume_rate_m3_s: f64,
        material: AdvectedMaterial,
    },
    /// Positive removed volume rate [m³/s]. The accepted local donor state
    /// supplies its density and h; extraction cannot prescribe a second h.
    Outflow { volume_rate_m3_s: f64 },
}

/// Authoritative extensive inventory and explicit sources on the original mesh.
#[derive(Clone, Copy, Debug)]
pub struct MaterialAdvectionInput<'a> {
    /// One positive mass per fluid cell [kg]; exactly zero in solid cells.
    /// Retain returned masses between steps instead of reconstructing rho*V.
    pub mass_kg: &'a [f64],
    /// Original specific enthalpy per fluid cell [J/kg]; solid entries ignored.
    pub specific_enthalpy_j_kg: &'a [f64],
    /// One entry per cell; solid entries must be `None`.
    pub volume_sources: &'a [MaterialVolumeSource],
    /// Incoming/backflow state for x-, x+, y-, y+, z-, z+ boundaries.
    /// Every actually incoming boundary needs a state. Outflow uses local h.
    pub boundary_inflow: [Option<AdvectedMaterial>; 6],
}

/// Pressure and material-transport budgets for one operator-split substep.
#[derive(Clone, Copy, Debug)]
pub struct MaterialAdvectionControls {
    /// The shared physical duration and original pressure/work budgets.
    pub projection: DensityProjectionControls,
    /// Maximum voxels, including solids, before allocation.
    pub max_cells: usize,
    /// Maximum nonzero internal, boundary and source volume routes.
    pub max_connections: usize,
    /// Native Krylov iteration ceiling for EACH donor transport solve.
    pub max_iterations: usize,
    /// Relative true-residual tolerance in the original row-scaled solver.
    pub relative_tolerance: f64,
    /// Sum of absolute local inventory/donor-equation defects plus the global
    /// boundary/source mass defect [kg]. No mass repair is performed.
    pub maximum_mass_defect_kg: f64,
    /// Sum of absolute local mass*h equation defects plus the independently
    /// reduced global enthalpy defect [J]. No numerical defect becomes heat.
    pub maximum_enthalpy_defect_j: f64,
}

/// Readback of actual represented inventories and their physical transfers.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialAdvectionReport {
    pub mass_before_kg: f64,
    pub mass_after_kg: f64,
    /// Boundary inflow plus local injected material [kg].
    pub incoming_mass_kg: f64,
    /// Boundary outflow plus local extracted material [kg].
    pub outgoing_mass_kg: f64,
    /// Signed enthalpy carried across those boundaries/sources [J].
    pub carried_enthalpy_j: f64,
    /// Actual inventory enthalpy change, reduced around one datum [J].
    pub enthalpy_change_j: f64,
    pub absolute_mass_defect_kg: f64,
    pub absolute_enthalpy_defect_j: f64,
    pub density_iterations: usize,
    pub enthalpy_iterations: usize,
    pub mass_routes: usize,
}

/// One complete pressure/material-transport candidate, without input mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialAdvectionStep {
    pub projection: DensityProjection,
    /// Conserved endpoint inventories; density is mass/domain volume.
    pub mass_kg: Vec<f64>,
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub report: MaterialAdvectionReport,
    mass_transport: AcceptedMassTransport,
}

impl MaterialAdvectionStep {
    /// Original accepted transfers, never reconstructed from a returned density.
    /// The endpoint inventory is public for existing consumers; reject a caller
    /// edit before using its private, immutable transport operator.
    pub(super) fn accepted_mass_transport(
        &self,
        gate: &CancelGate,
    ) -> Result<&AcceptedMassTransport, ChtError> {
        let transport = &self.mass_transport;
        if self.mass_kg.len() != transport.row_of.len() {
            return Err(refuse(
                "accepted inventory shape was changed before dependent transport",
            ));
        }
        for (c, &mass) in self.mass_kg.iter().enumerate() {
            poll_every(gate, c)?;
            let row = transport.row_of[c];
            let expected = if row == usize::MAX {
                0.0
            } else {
                transport.endpoint_mass_kg[row]
            };
            if mass != expected {
                return Err(refuse(
                    "accepted mass inventory was changed before dependent transport",
                ));
            }
        }
        Ok(transport)
    }
}

impl FlowField {
    /// Project volume continuity and advect actual mass and mass*h on the SAME
    /// native face fluxes. Positive volume sources inject their declared rho/h;
    /// extraction and boundary outflow carry the accepted implicit donor state.
    ///
    /// Density uses backward-Euler upwinding with cell volumes as storage.
    /// Its solved donor densities define transferred masses once. Paired
    /// transfers produce the authoritative endpoint masses, which are read
    /// back against the original density equation. A second native sparse
    /// solve carries mass*h over exactly these transfers, centered separately
    /// in each connected route component to exclude an arbitrary caloric datum
    /// from the solver tolerance. Actual stored h passes an independent audit.
    ///
    /// This is an operator-split pressure/material-transport component for a
    /// single-velocity density field. The caller supplies the momentum
    /// predictor, source momentum/work and a compatible material/EOS closure.
    /// It does not turn a thermodynamic density change into a material source,
    /// or infer phase kinetics, thermal expansion, mixing heat or pressure heat.
    /// First-order donor diffusion and temporal splitting require refinement.
    /// Pressure correction loss stays in the original mechanical account.
    /// Domain identity has the original projection's caller obligation.
    ///
    /// All budgets, positive inventories, actual mass/enthalpy equations and
    /// cancellation must pass before either corrected flow or material state
    /// is returned. Negative mass is never clipped or redistributed.
    pub fn project_and_advect_material(
        &self,
        domain: &VoxelDomain,
        faces: [FlowFace; 6],
        input: MaterialAdvectionInput<'_>,
        controls: MaterialAdvectionControls,
        gate: &CancelGate,
    ) -> Result<MaterialAdvectionStep, ChtError> {
        let density = admit(domain, input, controls, gate)?;
        let source: Vec<_> = input
            .volume_sources
            .iter()
            .map(|s| match s {
                MaterialVolumeSource::None => 0.0,
                MaterialVolumeSource::Inflow {
                    volume_rate_m3_s, ..
                } => *volume_rate_m3_s,
                MaterialVolumeSource::Outflow { volume_rate_m3_s } => -*volume_rate_m3_s,
            })
            .collect();
        let projected = self.project_density_volume_source(
            domain,
            faces,
            &density,
            &source,
            controls.projection,
            gate,
        )?;
        advect_projected_material(domain, projected, input, controls, gate)
    }
}

fn material(state: AdvectedMaterial) -> Result<(), ChtError> {
    finite_positive("flow.incoming_density", state.density_kg_m3)?;
    finite("flow.incoming_enthalpy", state.specific_enthalpy_j_kg)
}

pub(super) fn admit(
    domain: &VoxelDomain,
    input: MaterialAdvectionInput<'_>,
    controls: MaterialAdvectionControls,
    gate: &CancelGate,
) -> Result<Vec<f64>, ChtError> {
    poll(gate)?;
    let n = domain.cell_count();
    if n > controls.max_cells
        || input.mass_kg.len() != n
        || input.specific_enthalpy_j_kg.len() != n
        || input.volume_sources.len() != n
        || controls.max_iterations == 0
    {
        return Err(refuse(
            "inventory shape or transport work budget is invalid",
        ));
    }
    finite_positive("flow.material_duration", controls.projection.duration_s)?;
    finite_positive("flow.material_tolerance", controls.relative_tolerance)?;
    for bound in [
        controls.maximum_mass_defect_kg,
        controls.maximum_enthalpy_defect_j,
    ] {
        if !bound.is_finite() || bound < 0.0 {
            return Err(refuse("invalid material balance budget"));
        }
    }
    for value in input.boundary_inflow.into_iter().flatten() {
        material(value)?;
    }
    let mut density = vec![0.0; n];
    for c in 0..n {
        poll_every(gate, c)?;
        if !domain.is_fluid(c) {
            if input.mass_kg[c] != 0.0 || input.volume_sources[c] != MaterialVolumeSource::None {
                return Err(refuse(
                    "solid cells cannot contain advected inventory or sources",
                ));
            }
            continue;
        }
        finite_positive("flow.cell_mass", input.mass_kg[c])?;
        finite("flow.cell_enthalpy", input.specific_enthalpy_j_kg[c])?;
        density[c] = input.mass_kg[c] / domain.volume(c);
        finite_positive("flow.inventory_density", density[c])?;
        match input.volume_sources[c] {
            MaterialVolumeSource::None => {}
            MaterialVolumeSource::Inflow {
                volume_rate_m3_s,
                material: value,
            } => {
                finite_positive("flow.injected_volume_rate", volume_rate_m3_s)?;
                material(value)?;
            }
            MaterialVolumeSource::Outflow { volume_rate_m3_s } => {
                finite_positive("flow.removed_volume_rate", volume_rate_m3_s)?;
            }
        }
    }
    Ok(density)
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum InletLocation {
    Boundary(usize),
    Source(usize),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Route {
    pub(super) from: Option<usize>,
    pub(super) to: Option<usize>,
    volume: f64,
    pub(super) inlet: Option<AdvectedMaterial>,
    inlet_location: Option<InletLocation>,
    pub(super) mass: f64,
}

/// The same accepted mass routes and donor operator serve momentum and EOS
/// consumers. Private construction prevents an independently authored flux from
/// being passed off as the mass transport that produced the returned inventory.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct AcceptedMassTransport {
    pub(super) cells: Vec<usize>,
    row_of: Vec<usize>,
    pub(super) initial_mass_kg: Vec<f64>,
    pub(super) endpoint_mass_kg: Vec<f64>,
    pub(super) routes: Vec<Route>,
    parents: Vec<usize>,
    closure: Vec<f64>,
    operator: Csr,
    diagonal: Vec<f64>,
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    compensation: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let corrected = value - self.compensation;
        let next = self.value + corrected;
        self.compensation = (next - self.value) - corrected;
        self.value = next;
    }
}

fn root(parents: &mut [usize], mut i: usize) -> usize {
    while parents[i] != i {
        parents[i] = parents[parents[i]];
        i = parents[i];
    }
    i
}

/// Numerical transport on an already projected field. Keeping this separate
/// lets an admitted thermodynamic owner supply a volume change without falsely
/// declaring material injection. The enclosing owner owns that volume closure.
pub(super) fn advect_projected_material(
    domain: &VoxelDomain,
    projected: DensityProjection,
    input: MaterialAdvectionInput<'_>,
    controls: MaterialAdvectionControls,
    gate: &CancelGate,
) -> Result<MaterialAdvectionStep, ChtError> {
    let initial_density = admit(domain, input, controls, gate)?;
    let cells: Vec<_> = (0..domain.cell_count())
        .filter(|&c| domain.is_fluid(c))
        .collect();
    if cells.is_empty() {
        return Err(refuse("material advection needs fluid inventory"));
    }
    let n = cells.len();
    let mut row_of = vec![usize::MAX; domain.cell_count()];
    for (row, &c) in cells.iter().enumerate() {
        row_of[c] = row;
    }
    let dt = controls.projection.duration_s;
    let mut routes = Vec::new();
    let mut push = |from, to, volume: f64, inlet, inlet_location| -> Result<(), ChtError> {
        finite_positive("flow.transferred_volume", volume)?;
        if routes.len() >= controls.max_connections {
            return Err(refuse("material connection budget exhausted"));
        }
        routes.push(Route {
            from,
            to,
            volume,
            inlet,
            inlet_location,
            mass: 0.0,
        });
        Ok(())
    };
    for (row, &c) in cells.iter().enumerate() {
        poll_every(gate, row)?;
        for f in 0..6 {
            let nb = domain.neighbor(c, f);
            if nb.is_some_and(|other| f % 2 == 0 || !domain.is_fluid(other)) {
                continue;
            }
            let q = projected.field.outward(domain, c, f);
            finite("flow.material_face_flux", q)?;
            if q == 0.0 {
                continue;
            }
            let other = nb.map(|other| row_of[other]);
            if q > 0.0 {
                push(Some(row), other, dt * q, None, None)?;
            } else {
                let inlet =
                    if nb.is_none() {
                        Some(input.boundary_inflow[f].ok_or_else(|| {
                            refuse("incoming boundary lacks density and enthalpy")
                        })?)
                    } else {
                        None
                    };
                push(
                    other,
                    Some(row),
                    dt * -q,
                    inlet,
                    nb.is_none().then_some(InletLocation::Boundary(f)),
                )?;
            }
        }
        match input.volume_sources[c] {
            MaterialVolumeSource::None => {}
            MaterialVolumeSource::Inflow {
                volume_rate_m3_s,
                material,
            } => {
                push(
                    None,
                    Some(row),
                    dt * volume_rate_m3_s,
                    Some(material),
                    Some(InletLocation::Source(c)),
                )?;
            }
            MaterialVolumeSource::Outflow { volume_rate_m3_s } => {
                push(Some(row), None, dt * volume_rate_m3_s, None, None)?;
            }
        }
    }
    let mut diagonal: Vec<_> = cells.iter().map(|&c| domain.volume(c)).collect();
    let mut rhs: Vec<_> = cells.iter().map(|&c| input.mass_kg[c]).collect();
    let mut rho: Vec<_> = cells.iter().map(|&c| initial_density[c]).collect();
    let mut coo = Coo::new(n, n);
    let mut parents: Vec<_> = (0..n).collect();
    for (index, route) in routes.iter().enumerate() {
        poll_every(gate, index)?;
        if let Some(from) = route.from {
            diagonal[from] += route.volume;
        }
        match (route.from, route.to) {
            (Some(from), Some(to)) => {
                coo.push(to, from, -route.volume);
                let a = root(&mut parents, from);
                let b = root(&mut parents, to);
                if a != b {
                    parents[b] = a;
                }
            }
            (None, Some(to)) => {
                rhs[to] += route.volume * route.inlet.expect("admitted inlet").density_kg_m3
            }
            _ => {}
        }
    }
    let (a, density_iterations) = solve(coo, &diagonal, &mut rhs, &mut rho, controls, gate)?;
    let mut deltas = vec![Sum::default(); n];
    let mut incoming = Sum::default();
    let mut outgoing = Sum::default();
    for (index, route) in routes.iter_mut().enumerate() {
        poll_every(gate, index)?;
        let donor = route.from.map_or_else(
            || route.inlet.expect("admitted inlet").density_kg_m3,
            |from| rho[from],
        );
        finite_positive("flow.accepted_donor_density", donor)?;
        route.mass = route.volume * donor;
        finite_positive("flow.transferred_mass", route.mass)?;
        if let Some(from) = route.from {
            deltas[from].add(-route.mass);
        } else {
            incoming.add(route.mass);
        }
        if let Some(to) = route.to {
            deltas[to].add(route.mass);
        } else {
            outgoing.add(route.mass);
        }
    }
    let mut mass = vec![0.0; domain.cell_count()];
    let mut closure = vec![0.0; n];
    let mut mass_defect = Sum::default();
    let mut total_before = Sum::default();
    let mut total_after = Sum::default();
    let mut actual_mass_change = Sum::default();
    for (row, &c) in cells.iter().enumerate() {
        poll_every(gate, row)?;
        mass[c] = input.mass_kg[c] + deltas[row].value;
        finite_positive("flow.accepted_cell_mass", mass[c])?;
        rho[row] = mass[c] / domain.volume(c);
        finite_positive("flow.accepted_density", rho[row])?;
        let actual = mass[c] - input.mass_kg[c];
        closure[row] = deltas[row].value - actual;
        mass_defect.add(closure[row].abs());
        actual_mass_change.add(actual);
        total_before.add(input.mass_kg[c]);
        total_after.add(mass[c]);
    }
    mass_defect.add(residual(&a, &rhs, &rho, &diagonal, gate)?);
    let net_external_mass = incoming.value - outgoing.value;
    let global_mass_defect = actual_mass_change.value - net_external_mass;
    mass_defect.add(global_mass_defect.abs());
    finite("flow.material_mass_defect", mass_defect.value)?;
    if mass_defect.value > controls.maximum_mass_defect_kg {
        return Err(refuse(
            "actual mass transport exceeds its independent balance budget",
        ));
    }

    let mut references = vec![None; n];
    for row in 0..n {
        poll_every(gate, row)?;
        parents[row] = root(&mut parents, row);
        references[parents[row]].get_or_insert(input.specific_enthalpy_j_kg[cells[row]]);
    }
    let reference: Vec<_> = (0..n)
        .map(|row| references[parents[row]].expect("positive old inventory"))
        .collect();
    let mut h: Vec<_> = cells
        .iter()
        .enumerate()
        .map(|(row, &c)| input.specific_enthalpy_j_kg[c] - reference[row])
        .collect();
    let mut rhs: Vec<_> = cells
        .iter()
        .enumerate()
        .map(|(row, &c)| closure[row].mul_add(reference[row], input.mass_kg[c] * h[row]))
        .collect();
    let mut diagonal: Vec<_> = cells.iter().map(|&c| mass[c]).collect();
    let mut coo = Coo::new(n, n);
    for (index, route) in routes.iter().enumerate() {
        poll_every(gate, index)?;
        if let Some(from) = route.from {
            diagonal[from] += route.mass;
        }
        match (route.from, route.to) {
            (Some(from), Some(to)) => coo.push(to, from, -route.mass),
            (None, Some(to)) => {
                rhs[to] += route.mass
                    * (route.inlet.expect("admitted inlet").specific_enthalpy_j_kg - reference[to])
            }
            _ => {}
        }
    }
    let (a, enthalpy_iterations) = solve(coo, &diagonal, &mut rhs, &mut h, controls, gate)?;
    let mut endpoint_h = vec![0.0; domain.cell_count()];
    for (row, &c) in cells.iter().enumerate() {
        poll_every(gate, row)?;
        endpoint_h[c] = h[row] + reference[row];
        finite("flow.accepted_enthalpy", endpoint_h[c])?;
        h[row] = endpoint_h[c] - reference[row];
    }
    let mut energy_defect = residual(&a, &rhs, &h, &diagonal, gate)?;
    let datum = input.specific_enthalpy_j_kg[cells[0]];
    let mut change = Sum::default();
    let mut carried = Sum::default();
    for (row, &c) in cells.iter().enumerate() {
        poll_every(gate, row)?;
        change.add(mass[c] * (endpoint_h[c] - input.specific_enthalpy_j_kg[c]));
        change.add((mass[c] - input.mass_kg[c]) * (input.specific_enthalpy_j_kg[c] - datum));
    }
    for (index, route) in routes.iter().enumerate() {
        poll_every(gate, index)?;
        if route.from.is_none() {
            carried.add(
                route.mass * (route.inlet.expect("admitted inlet").specific_enthalpy_j_kg - datum),
            );
        } else if route.to.is_none() {
            carried
                .add(-route.mass * (endpoint_h[cells[route.from.expect("outflow donor")]] - datum));
        }
    }
    energy_defect += (global_mass_defect.mul_add(datum, change.value - carried.value)).abs();
    let enthalpy_change_j = actual_mass_change.value.mul_add(datum, change.value);
    let carried_enthalpy_j = net_external_mass.mul_add(datum, carried.value);
    for value in [
        energy_defect,
        enthalpy_change_j,
        carried_enthalpy_j,
        total_before.value,
        total_after.value,
        incoming.value,
        outgoing.value,
    ] {
        finite("flow.material_account", value)?;
    }
    if energy_defect > controls.maximum_enthalpy_defect_j {
        return Err(refuse(
            "actual enthalpy transport exceeds its independent balance budget",
        ));
    }
    poll(gate)?;
    let mass_routes = routes.len();
    let mass_transport = AcceptedMassTransport {
        initial_mass_kg: cells.iter().map(|&c| input.mass_kg[c]).collect(),
        endpoint_mass_kg: cells.iter().map(|&c| mass[c]).collect(),
        cells,
        row_of,
        routes,
        parents,
        closure,
        operator: a,
        diagonal,
    };
    Ok(MaterialAdvectionStep {
        projection: projected,
        mass_kg: mass,
        specific_enthalpy_j_kg: endpoint_h,
        mass_transport,
        report: MaterialAdvectionReport {
            mass_before_kg: total_before.value,
            mass_after_kg: total_after.value,
            incoming_mass_kg: incoming.value,
            outgoing_mass_kg: outgoing.value,
            carried_enthalpy_j,
            enthalpy_change_j,
            absolute_mass_defect_kg: mass_defect.value,
            absolute_enthalpy_defect_j: energy_defect,
            density_iterations,
            enthalpy_iterations,
            mass_routes,
        },
    })
}

fn solve(
    mut coo: Coo,
    diagonal: &[f64],
    rhs: &mut [f64],
    x: &mut [f64],
    controls: MaterialAdvectionControls,
    gate: &CancelGate,
) -> Result<(Csr, usize), ChtError> {
    for row in 0..diagonal.len() {
        poll_every(gate, row)?;
        finite_positive("flow.material_storage", diagonal[row])?;
        finite("flow.material_rhs", rhs[row])?;
        finite("flow.material_guess", x[row])?;
        coo.push(row, row, diagonal[row]);
    }
    let a = scale_rows(&coo, rhs);
    let result = bicgstab_ilu0(
        "material donor transport",
        &a,
        rhs,
        x,
        controls.relative_tolerance,
        controls.max_iterations,
        gate,
    )?;
    Ok((a, result.iterations))
}

fn residual(
    a: &Csr,
    rhs: &[f64],
    x: &[f64],
    diagonal: &[f64],
    gate: &CancelGate,
) -> Result<f64, ChtError> {
    let mut defect = Sum::default();
    for row in 0..x.len() {
        poll_every(gate, row)?;
        let (columns, values) = a.row(row);
        let mut sum = Sum::default();
        sum.add(-rhs[row]);
        for (&column, &value) in columns.iter().zip(values) {
            sum.add(value * x[column]);
        }
        defect.add((diagonal[row] * sum.value).abs());
    }
    finite("flow.material_readback", defect.value)?;
    Ok(defect.value)
}

fn refuse(reason: &str) -> ChtError {
    ChtError::InvalidInput {
        field: "flow.material_advection",
        reason: reason.into(),
    }
}

#[cfg(test)]
mod tests;
