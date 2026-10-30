//! Electroquasistatic charge relaxation on the original dielectric/Ohmic mesh.
//!
//! `(K_epsilon1 V1-K_epsilon0 V0)/dt + K_sigma (V1+V0)/2 = I_external` uses the existing
//! material FEEC operators and checked scalar-potential solve. Local midpoint
//! E.J comes from the original Ohmic recovery; endpoint E.D/2 comes from the
//! original dielectric recovery. No fitted RC element or duplicate cell law is
//! introduced. Native nodal charges include volume and boundary weak charges.
//!
//! Ordinary construction requires source-declared constant permittivity. An
//! explicitly bound thermal derivative law also admits temperature-dependent
//! storage, with local signed parameter work exposed for its caloric owner.
//! Conductivity and epsilon retain the same complete material bundle. Geometry and voltage
//! trace topology remain fixed. Unprescribed faces carry zero TOTAL current,
//! including displacement current. This is a quasistatic leaky dielectric, not
//! dielectric dispersion, magnetic induction, breakdown or carrier transport.

use super::{ElectrostaticError, PlanarElectrostatics};
use crate::conduction::{ConductionError, PlanarConduction};
use crate::scalar_potential::{self, PotentialError};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_material::{
    conductor::{ResolvedConductorState, planar::PlanarConductorMaterial},
    field_storage::{planar::PlanarDielectricLaw,
        thermal::electric::ThermalPlanarDielectricMaterial},
};
use fs_solver::{CheckedCgConfig, SolveReport};
use fs_sparse::Csr;
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
pub struct DielectricRelaxationControls {
    pub linear: CheckedCgConfig,
    pub charge_tolerance_c: f64,
    pub current_tolerance_a: f64,
    pub energy_tolerance_j: f64,
}

#[derive(Debug)]
pub enum DielectricRelaxationError {
    Electric(ElectrostaticError),
    Conduction(ConductionError),
    Invalid(&'static str),
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for DielectricRelaxationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native dielectric charge relaxation: {self:?}")
    }
}
impl std::error::Error for DielectricRelaxationError {}
impl From<ElectrostaticError> for DielectricRelaxationError {
    fn from(value: ElectrostaticError) -> Self {
        match value {
            ElectrostaticError::Cancelled => Self::Cancelled,
            other => Self::Electric(other),
        }
    }
}
impl From<ConductionError> for DielectricRelaxationError {
    fn from(value: ConductionError) -> Self {
        match value {
            ConductionError::Cancelled => Self::Cancelled,
            other => Self::Conduction(other),
        }
    }
}
impl From<PotentialError> for DielectricRelaxationError {
    fn from(value: PotentialError) -> Self {
        ElectrostaticError::from(value).into()
    }
}
type Error = DielectricRelaxationError;
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid("dielectric arithmetic is unrepresentable"))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_controls(controls: DielectricRelaxationControls) -> Result<(), Error> {
    if !controls.linear.tolerance.is_finite()
        || controls.linear.tolerance <= 0.0
        || controls.linear.tolerance >= 1.0
        || [
            controls.charge_tolerance_c,
            controls.current_tolerance_a,
            controls.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err(Error::Invalid(
            "positive finite native solve and physical balance tolerances required",
        ));
    }
    Ok(())
}

/// Accepted electric potential/charge, fixed dielectric storage and actual
/// scalar or oriented conductivity. Each conductor retains the COMPLETE common
/// source bundle, including every selected absolute permittivity coefficient.
/// Dielectric material axes are selected explicitly and held independently of
/// the conductor axes; both resolve into the same named mesh-coordinate frame.
#[derive(Clone)]
pub struct PlanarLeakyDielectric<C = ResolvedConductorState> {
    dielectric: Arc<PlanarElectrostatics>,
    thermal_laws: Option<Arc<[ThermalPlanarDielectricMaterial]>>,
    conduction: Arc<PlanarConduction<C>>,
    voltage_v: Vec<f64>,
    nodal_charge_c: Vec<f64>,
    cell_energy_j: Vec<f64>,
    stored_energy_j: f64,
    voltage_vertices: Vec<usize>,
    time_s: f64,
}

/// Local Helmholtz storage and explicit coefficient work on one original cell.
/// P=-volume*E0.(epsilon1-epsilon0).E1/2 is the exact two-endpoint work required
/// by this native charge equation. It is reversible transfer, never Joule heat.
#[derive(Clone, Debug, PartialEq)]
pub struct DielectricEnergyCell {
    pub triangle: usize,
    pub displacement_c_m2: [[f64; 2]; 2],
    pub free_energy_j: [f64; 2],
    pub parameter_work_j: f64,
}

/// One physical midpoint interval. Current densities are resolved locally in
/// mesh coordinates; source currents contain conduction AND displacement.
#[derive(Clone, Debug)]
pub struct DielectricRelaxationStep {
    pub interval_s: [f64; 2],
    pub midpoint_voltage_v: Vec<f64>,
    pub electric_field_v_m: Vec<[f64; 2]>,
    pub conduction_current_density_a_m2: Vec<[f64; 2]>,
    pub displacement_current_density_a_m2: Vec<[f64; 2]>,
    pub cell_power_w: Vec<f64>,
    pub nodal_boundary_current_a: Vec<Option<f64>>,
    pub maximum_free_current_residual_a: f64,
    pub boundary_current_residual_a: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub heat_j: f64,
    pub supplied_work_j: f64,
    pub parameter_work_j: f64,
    pub electric_cells: Vec<DielectricEnergyCell>,
    pub energy_residual_j: f64,
    pub linear: Option<SolveReport>,
}

pub struct DielectricRelaxationCandidate<C = ResolvedConductorState> {
    field: PlanarLeakyDielectric<C>,
    step: DielectricRelaxationStep,
}
impl<C> DielectricRelaxationCandidate<C> {
    pub fn field(&self) -> &PlanarLeakyDielectric<C> {
        &self.field
    }
    pub fn step(&self) -> &DielectricRelaxationStep {
        &self.step
    }
    pub fn into_parts(self) -> (PlanarLeakyDielectric<C>, DielectricRelaxationStep) {
        (self.field, self.step)
    }
}

fn dielectric<C: PlanarConductorMaterial>(
    conduction: &PlanarConduction<C>,
    laws: Vec<PlanarDielectricLaw>,
    thermal_laws: Option<&[ThermalPlanarDielectricMaterial]>,
    limits: WeightedAssemblyLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PlanarElectrostatics, Error> {
    if !conduction.contacts().is_empty() {
        return Err(Error::Invalid(
            "electrical contact charge and heat evolution require an owned surface law",
        ));
    }
    if laws.len() != conduction.conductors().len() || laws.len() > limits.max_cells
        || thermal_laws.is_some_and(|thermal| thermal.len() != laws.len()) {
        return Err(Error::Invalid("one bounded dielectric source law per original conductor cell required"));
    }
    let mut geometry_frame = conduction.conductors().iter()
        .find_map(|material| material.frame().map(|frame| frame.geometry_frame()));
    let mut materials = Vec::with_capacity(laws.len());
    for (triangle, (conductor, law)) in conduction.conductors().iter().zip(&laws).enumerate() {
        poll(cancelled)?;
        if let Some(frame) = law.frame() {
            if geometry_frame.is_some_and(|id| id != frame.geometry_frame()) {
                return Err(Error::Invalid(
                    "dielectric and conductivity frames must refer to the same mesh coordinates",
                ));
            }
            geometry_frame = Some(frame.geometry_frame());
        }
        let source = conductor.material();
        if let Some(thermal) = thermal_laws {
            if thermal[triangle].law() != *law {
                return Err(Error::Invalid("dielectric derivative law changed original material axes"));
            }
            thermal[triangle].at_material(source).map_err(ElectrostaticError::Material)?;
        } else {
            law.resolve_constant(source).map_err(ElectrostaticError::Material)?;
        }
        materials.push(source.clone());
    }
    Ok(PlanarElectrostatics::new_with_laws(
        conduction.mesh().clone(),
        materials,
        laws,
        limits,
        cancelled,
    )?)
}

// Difference form preserves gauge invariance before multiplying by a Hodge.
fn apply(
    matrix: &Csr,
    voltage: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, Error> {
    let mut result = Vec::with_capacity(matrix.nrows());
    for row in 0..matrix.nrows() {
        poll(cancelled)?;
        let (columns, values) = matrix.row(row);
        let mut total = 0.0;
        for (&column, &value) in columns.iter().zip(values) {
            total = finite(value.mul_add(finite(voltage[column] - voltage[row])?, total))?;
        }
        result.push(total);
    }
    Ok(result)
}

// Merge the existing sparse operators; no second FE discretization is assembled.
fn midpoint_operator(
    epsilon: &Csr,
    sigma: &Csr,
    half_dt: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Csr, Error> {
    let n = epsilon.nrows();
    let mut row_ptr = vec![0];
    let mut columns = Vec::new();
    let mut values = Vec::new();
    for row in 0..n {
        poll(cancelled)?;
        let (ec, ev) = epsilon.row(row);
        let (sc, sv) = sigma.row(row);
        let (mut i, mut j) = (0, 0);
        while i < ec.len() || j < sc.len() {
            let column = ec
                .get(i)
                .copied()
                .unwrap_or(usize::MAX)
                .min(sc.get(j).copied().unwrap_or(usize::MAX));
            let e = if ec.get(i) == Some(&column) {
                let x = ev[i];
                i += 1;
                x
            } else {
                0.0
            };
            let s = if sc.get(j) == Some(&column) {
                let x = sv[j];
                j += 1;
                x
            } else {
                0.0
            };
            columns.push(column);
            values.push(finite(half_dt.mul_add(s, e))?);
        }
        row_ptr.push(columns.len());
    }
    Ok(Csr::from_parts(n, n, row_ptr, columns, values))
}

fn admit_refresh<C: PlanarConductorMaterial>(old: &[C], new: &[C]) -> Result<(), Error> {
    if old.len() != new.len() {
        return Err(Error::Invalid(
            "original dielectric/conductor coverage changed",
        ));
    }
    for (a, b) in old.iter().zip(new) {
        let (before, after) = (a.material(), b.material());
        if !a.has_same_electrical_source(b)
            || a.frame() != b.frame()
            || a.liquid_mass_fraction_axis() != b.liquid_mass_fraction_axis()
            || before.card_identity() != after.card_identity()
            || before.material() != after.material()
            || before.axis_quantities() != after.axis_quantities()
            || before
                .query_point()
                .iter()
                .map(|(name, _)| name)
                .ne(after.query_point().iter().map(|(name, _)| name))
            || before.properties().len() != after.properties().len()
            || before
                .properties()
                .iter()
                .zip(after.properties())
                .any(|(p, q)| {
                    p.requirement() != q.requirement()
                        || p.source_axes() != q.source_axes()
                        || p.answer().receipt.selected != q.answer().receipt.selected
                })
        {
            return Err(Error::Invalid(
                "dielectric refresh must retain original source bundles, basis, frames and coordinate meanings",
            ));
        }
    }
    Ok(())
}

impl<C: PlanarConductorMaterial> PlanarLeakyDielectric<C> {
    /// Initialize a physical P1 electric-potential state. Its weak free charges
    /// are K_epsilon V, including surface charge; they are not an additional
    /// prescribed source. `voltage_vertices` selects complete exterior traces
    /// held at their initial values. Other boundaries subsequently isolate total
    /// current. Every component needs an anchored trace. A precharged field may
    /// be obtained from the existing electrostatic solve without fitting a C.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        conduction: PlanarConduction<C>,
        voltage_v: Vec<f64>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: DielectricRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if conduction.conductors().len() > limits.max_cells {
            return Err(Error::Invalid("dielectric cell count exceeds assembly budget"));
        }
        let laws = vec![PlanarDielectricLaw::Isotropic; conduction.conductors().len()];
        Self::new_with_dielectric_laws(
            conduction, laws, voltage_v, voltage_vertices, time_s, limits, controls, cancelled,
        )
    }

    /// Initialize the SAME charge owner with explicit source dielectric laws.
    /// Each coefficient is selected from its conductor's complete material
    /// bundle; no second material state or conductivity-derived axes are used.
    /// Scalar and oriented cells may coexist. Both principal permittivities,
    /// their selected claims and their independent frames remain fixed through
    /// subsequent conductivity/temperature/phase refreshes.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_dielectric_laws(
        conduction: PlanarConduction<C>,
        laws: Vec<PlanarDielectricLaw>,
        voltage_v: Vec<f64>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: DielectricRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_with_storage(conduction, laws, None, voltage_v, voltage_vertices,
            time_s, limits, controls, cancelled)
    }

    /// Initialize the same charge owner with source-backed epsilon(T) and its
    /// exact thermal derivatives. Every triangle retains its original bundle
    /// and fixed axes. Subsequent native steps expose parameter work; a thermal
    /// consumer must additionally close field entropy and caloric internal
    /// energy, including positive combined capacity over its field envelope.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_thermal_dielectric_laws(
        conduction: PlanarConduction<C>,
        thermal_laws: Vec<ThermalPlanarDielectricMaterial>,
        voltage_v: Vec<f64>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: DielectricRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if thermal_laws.len() > limits.max_cells {
            return Err(Error::Invalid("dielectric thermal laws exceed assembly budget"));
        }
        let laws = thermal_laws.iter().map(ThermalPlanarDielectricMaterial::law).collect();
        Self::new_with_storage(conduction, laws, Some(Arc::from(thermal_laws)), voltage_v,
            voltage_vertices, time_s, limits, controls, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn new_with_storage(
        conduction: PlanarConduction<C>,
        laws: Vec<PlanarDielectricLaw>,
        thermal_laws: Option<Arc<[ThermalPlanarDielectricMaterial]>>,
        voltage_v: Vec<f64>,
        voltage_vertices: &[usize],
        time_s: f64,
        limits: WeightedAssemblyLimits,
        controls: DielectricRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        admit_controls(controls)?;
        if !time_s.is_finite()
            || time_s < 0.0
            || voltage_v.is_empty()
            || voltage_v.len() != conduction.mesh().vertices().len()
            || voltage_v.iter().any(|v| !v.is_finite())
            || voltage_vertices.len() > voltage_v.len()
            || voltage_vertices.iter().any(|v| *v >= voltage_v.len())
        {
            return Err(Error::Invalid(
                "finite time, complete potentials and bounded voltage traces required",
            ));
        }
        let boundary: Vec<_> = voltage_vertices
            .iter()
            .map(|&v| (v, voltage_v[v]))
            .collect();
        scalar_potential::prescribed_boundary(conduction.mesh(), &boundary, cancelled)?;
        let dielectric = dielectric(&conduction, laws, thermal_laws.as_deref(), limits, cancelled)?;
        let cells = dielectric.cell_response(&voltage_v, cancelled)?;
        let nodal_charge_c = apply(&dielectric.stiffness, &voltage_v, cancelled)?;
        check_storage(&voltage_v, &nodal_charge_c, cells.energy_j, controls)?;
        let mut voltage_vertices = voltage_vertices.to_vec();
        voltage_vertices.sort_unstable();
        poll(cancelled)?;
        Ok(Self {
            dielectric: Arc::new(dielectric),
            thermal_laws,
            conduction: Arc::new(conduction),
            voltage_v,
            nodal_charge_c,
            cell_energy_j: cells.cell_energy_j,
            stored_energy_j: cells.energy_j,
            voltage_vertices,
            time_s,
        })
    }
    pub fn conduction(&self) -> &PlanarConduction<C> {
        &self.conduction
    }
    pub fn dielectric(&self) -> &PlanarElectrostatics {
        &self.dielectric
    }
    pub fn thermal_dielectric_laws(&self) -> Option<&[ThermalPlanarDielectricMaterial]> {
        self.thermal_laws.as_deref()
    }
    /// Actual local D from the original native field recovery, in mesh axes.
    pub fn displacement_c_m2(
        &self, cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<[f64; 2]>, Error> {
        Ok(self.dielectric.cell_response(&self.voltage_v, cancelled)?.displacement_c_m2)
    }
    pub fn voltage_v(&self) -> &[f64] {
        &self.voltage_v
    }
    pub fn nodal_charge_c(&self) -> &[f64] {
        &self.nodal_charge_c
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Advance from the same accepted charge on every call. Endpoint voltage
    /// values describe a linear ramp over the interval on the original trace;
    /// the recovered total source current pairs with its midpoint potential.
    /// Refreshed conductivity is held over this interval. Standalone callers
    /// own its state history; thermal coupling must solve that same endpoint.
    /// Midpoint resolves held-coefficient dynamics to second order but does not
    /// preserve charge sign at arbitrary time steps; resolve the physical
    /// relaxation time and refine time/mesh independently of algebraic residuals.
    #[allow(clippy::too_many_arguments)]
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        conductors: Vec<C>,
        endpoint_voltage_v: &[(usize, f64)],
        limits: WeightedAssemblyLimits,
        controls: DielectricRelaxationControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DielectricRelaxationCandidate<C>, Error> {
        poll(cancelled)?;
        admit_controls(controls)?;
        let dt = finite(target_s - self.time_s)?;
        if dt <= 0.0 {
            return Err(Error::Invalid("positive dielectric interval required"));
        }
        admit_refresh(self.conduction.conductors(), &conductors)?;
        let conduction = PlanarConduction::new(
            self.conduction.mesh().clone(),
            conductors,
            limits,
            cancelled,
        )?;
        let dielectric = dielectric(
            &conduction, self.dielectric.permittivity_laws().to_vec(),
            self.thermal_laws.as_deref(), limits, cancelled,
        )?;
        if self.thermal_laws.is_none() && dielectric.permittivity_coefficients().iter()
            .zip(self.dielectric.permittivity_coefficients())
            .any(|(a, b)| a.matrix_f_m().into_iter().flatten()
                .zip(b.matrix_f_m().into_iter().flatten())
                .any(|(left, right)| left.to_bits() != right.to_bits()))
        {
            return Err(Error::Invalid(
                "dielectric permittivity changed without its parameter-work owner",
            ));
        }
        let prescribed = scalar_potential::prescribed_boundary(
            conduction.mesh(),
            endpoint_voltage_v,
            cancelled,
        )?;
        if prescribed
            .iter()
            .enumerate()
            .filter_map(|(v, x)| x.map(|_| v))
            .ne(self.voltage_vertices.iter().copied())
        {
            return Err(Error::Invalid(
                "live dielectric voltage-trace topology must remain fixed",
            ));
        }
        let half_dt = finite(0.5 * dt)?;
        if half_dt == 0.0 {
            return Err(Error::Invalid("dielectric interval underflow"));
        }
        let matrix = midpoint_operator(
            &dielectric.stiffness,
            conduction.stiffness_operator(),
            half_dt,
            cancelled,
        )?;
        let old_current = apply(conduction.stiffness_operator(), &self.voltage_v, cancelled)?;
        let rhs = self
            .nodal_charge_c
            .iter()
            .zip(old_current)
            .map(|(&q, i)| finite((-half_dt).mul_add(i, q)))
            .collect::<Result<Vec<_>, _>>()?;
        let (voltage_v, linear) =
            scalar_potential::solve(&matrix, &rhs, &prescribed, controls.linear, cancelled)?;
        let midpoint_voltage_v = self
            .voltage_v
            .iter()
            .zip(&voltage_v)
            .map(|(a, b)| finite(0.5 * a + 0.5 * b))
            .collect::<Result<Vec<_>, _>>()?;
        let losses = conduction.cell_response(&midpoint_voltage_v, cancelled)?;
        let previous = self.dielectric.cell_response(&self.voltage_v, cancelled)?;
        let endpoint = dielectric.cell_response(&voltage_v, cancelled)?;
        let mut electric_cells = Vec::with_capacity(previous.cell_energy_j.len());
        let mut parameter_work_j = 0.0;
        for triangle in 0..previous.cell_energy_j.len() {
            poll(cancelled)?;
            let before = self.dielectric.permittivity_coefficients()[triangle].matrix_f_m();
            let after = dielectric.permittivity_coefficients()[triangle].matrix_f_m();
            let e0 = previous.electric_field_v_m[triangle];
            let e1 = endpoint.electric_field_v_m[triangle];
            let mut delta_d = [0.0; 2];
            for i in 0..2 {
                for j in 0..2 {
                    delta_d[i] = finite(finite(after[i][j] - before[i][j])?.mul_add(e1[j], delta_d[i]))?;
                }
            }
            let volume = conduction.mesh().face_measure(triangle)
                .ok_or(Error::Invalid("dielectric cell lacks its original volume"))?;
            let work = finite(-0.5 * volume * finite(e0[0].mul_add(delta_d[0], e0[1] * delta_d[1]))?)?;
            parameter_work_j = finite(parameter_work_j + work)?;
            electric_cells.push(DielectricEnergyCell {
                triangle,
                displacement_c_m2: [previous.displacement_c_m2[triangle], endpoint.displacement_c_m2[triangle]],
                free_energy_j: [previous.cell_energy_j[triangle], endpoint.cell_energy_j[triangle]],
                parameter_work_j: work,
            });
        }
        let nodal_charge_c = apply(&dielectric.stiffness, &voltage_v, cancelled)?;
        check_storage(&voltage_v, &nodal_charge_c, endpoint.energy_j, controls)?;
        let current = apply(
            conduction.stiffness_operator(),
            &midpoint_voltage_v,
            cancelled,
        )?;
        let mut nodal_boundary_current_a = vec![None; voltage_v.len()];
        let mut maximum_free_current_residual_a: f64 = 0.0;
        let mut boundary_current_residual_a = 0.0;
        let mut source_power = 0.0;
        let reference = midpoint_voltage_v[self.voltage_vertices[0]];
        for (v, i) in current.into_iter().enumerate() {
            poll(cancelled)?;
            let total = finite(finite(nodal_charge_c[v] - self.nodal_charge_c[v])? / dt + i)?;
            if prescribed[v].is_some() {
                nodal_boundary_current_a[v] = Some(total);
                boundary_current_residual_a = finite(boundary_current_residual_a + total)?;
                source_power =
                    finite((midpoint_voltage_v[v] - reference).mul_add(total, source_power))?;
            } else {
                maximum_free_current_residual_a = maximum_free_current_residual_a.max(total.abs());
            }
        }
        let displacement_current_density_a_m2 = previous
            .displacement_c_m2
            .iter()
            .zip(&endpoint.displacement_c_m2)
            .map(|(a, b)| Ok([finite((b[0] - a[0]) / dt)?, finite((b[1] - a[1]) / dt)?]))
            .collect::<Result<Vec<_>, Error>>()?;
        let heat_j = finite(dt * losses.power_w)?;
        let supplied_work_j = finite(dt * source_power)?;
        let energy_residual_j =
            finite(endpoint.energy_j - self.stored_energy_j + heat_j - supplied_work_j - parameter_work_j)?;
        if (losses.power_w != 0.0 && heat_j == 0.0)
            || (source_power != 0.0 && supplied_work_j == 0.0)
            || maximum_free_current_residual_a > controls.current_tolerance_a
            || boundary_current_residual_a.abs() > controls.current_tolerance_a
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(Error::Balance(
                "charge continuity, source work and physical field/Joule energy do not close",
            ));
        }
        let step = DielectricRelaxationStep {
            interval_s: [self.time_s, target_s],
            midpoint_voltage_v,
            electric_field_v_m: losses.electric_field_v_m,
            conduction_current_density_a_m2: losses.current_density_a_m2,
            displacement_current_density_a_m2,
            cell_power_w: losses.cell_power_w,
            nodal_boundary_current_a,
            maximum_free_current_residual_a,
            boundary_current_residual_a,
            stored_before_j: self.stored_energy_j,
            stored_after_j: endpoint.energy_j,
            heat_j,
            supplied_work_j,
            parameter_work_j,
            electric_cells,
            energy_residual_j,
            linear,
        };
        poll(cancelled)?;
        Ok(DielectricRelaxationCandidate {
            field: Self {
                dielectric: Arc::new(dielectric),
                thermal_laws: self.thermal_laws.clone(),
                conduction: Arc::new(conduction),
                voltage_v,
                nodal_charge_c,
                cell_energy_j: endpoint.cell_energy_j,
                stored_energy_j: endpoint.energy_j,
                voltage_vertices: self.voltage_vertices.clone(),
                time_s: target_s,
            },
            step,
        })
    }
}

fn check_storage(
    voltage: &[f64],
    charge: &[f64],
    energy: f64,
    controls: DielectricRelaxationControls,
) -> Result<(), Error> {
    let reference = voltage[0];
    let mut pairing = 0.0;
    let mut total_charge = 0.0;
    for (&v, &q) in voltage.iter().zip(charge) {
        pairing = finite((v - reference).mul_add(q, pairing))?;
        total_charge = finite(total_charge + q)?;
    }
    if finite(energy - 0.5 * pairing)?.abs() > controls.energy_tolerance_j
        || total_charge.abs() > controls.charge_tolerance_c
    {
        return Err(Error::Balance(
            "original dielectric charge and cell energy disagree",
        ));
    }
    Ok(())
}
