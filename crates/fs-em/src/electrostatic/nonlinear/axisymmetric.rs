//! Reversible nonlinear electrostatics on the original revolved P1 complex.
//!
//! The shared Gauss/Newton owner integrates the actual meridional measure
//! `angle * r * dr * dz`. Electric fields are `[E_r, E_z]`; the azimuthal
//! component is zero. Constant cell fields have exact radial energy and
//! tangent integrals, while volume charge uses the original exact P1 radial
//! moments. Charge, energy and capacitance describe the authored angular span.
//!
//! Radial source laws and explicitly framed meridional tensors are admitted.
//! A fixed Cartesian anisotropy is not implicitly revolved. The angular span
//! is a symmetry sector, not independently driven azimuthal faces. There is
//! no inferred leakage, polarization history, breakdown or electrostriction.

use super::{
    DielectricMaterial, ElectrodeControl, Error, NonlinearElectricField, NonlinearElectrodes,
    NonlinearElectrostaticControls, NonlinearElectrostaticReport, PlanarNonlinearElectrostatics,
    add, finite, invalid, poll,
    tangent::{DifferentialCapacitance, DifferentialCapacitanceControls},
};
use crate::electrostatic::{
    ElectrostaticChargeSheet, ElectrostaticDisplacementTrace, PlanarElectrostaticSources,
    axisymmetric::AxisymmetricElectrostaticSources,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_rep_mesh::{Metric2, TriComplex2};

/// Original radial geometry, material sources and prepared matrix-free cells.
/// Its planar implementation detail cannot escape through planar force APIs.
pub struct AxisymmetricNonlinearElectrostatics {
    core: PlanarNonlinearElectrostatics,
}

impl AxisymmetricNonlinearElectrostatics {
    /// Prepare fixed-source linear/tensor and reversible radial D(E) cells.
    /// Tensor axes are in local `(r,z)` coordinates, with the omitted azimuthal
    /// direction decoupled. Every field trial retains original source support.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<DielectricMaterial>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Axisymmetric { .. }) {
            return Err(invalid(
                "axisymmetric nonlinear electrostatics requires an explicit radial metric",
            ));
        }
        #[cfg(feature = "dielectric-hysteresis")]
        if materials
            .iter()
            .any(|material| matches!(material, DielectricMaterial::Hysteretic(_)))
        {
            return Err(invalid(
                "axisymmetric retained polarization requires its own cylindrical history owner",
            ));
        }
        Ok(Self {
            core: PlanarNonlinearElectrostatics::assemble_metric(
                mesh, materials, limits, cancelled,
            )?,
        })
    }

    /// Original metre-valued `[radius, axial_coordinate]` complex.
    pub fn mesh(&self) -> &TriComplex2 {
        self.core.mesh()
    }

    /// Original complete source laws in cell order.
    pub fn materials(&self) -> &[DielectricMaterial] {
        self.core.materials()
    }

    /// Prepare disjoint physical electrodes, each a union of positive-area
    /// exterior surfaces of revolution. A symmetry-axis line is not an
    /// electrode, including when its requested total free charge is zero.
    pub fn prepare_electrodes<'m>(
        &'m self,
        terminals: &[&[usize]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricNonlinearElectrodes<'m>, Error> {
        Ok(AxisymmetricNonlinearElectrodes {
            problem: self,
            core: self.core.prepare_electrodes(terminals, cancelled)?,
        })
    }
}

/// Prepared equipotential quotient with original revolved surface identities.
pub struct AxisymmetricNonlinearElectrodes<'m> {
    problem: &'m AxisymmetricNonlinearElectrostatics,
    core: NonlinearElectrodes<'m>,
}

impl<'m> AxisymmetricNonlinearElectrodes<'m> {
    pub fn problem(&self) -> &'m AxisymmetricNonlinearElectrostatics {
        self.problem
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        self.core.terminals()
    }

    /// Solve div(D)=rho with voltage or TOTAL free-charge controls. A voltage
    /// reference is required on each electrically connected component. The
    /// initial potentials follow the shared component-relative convention.
    /// Natural unlisted surfaces have D.n=0; r=0 retains regularity.
    pub fn solve(
        &self,
        free_charge_density_c_m3: &[f64],
        controls: &[ElectrodeControl],
        initial_relative_voltage_v: &[f64],
        solve: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricNonlinearElectricField<'m>, Error> {
        Ok(AxisymmetricNonlinearElectricField {
            problem: self.problem,
            field: self.core.solve(
                free_charge_density_c_m3,
                controls,
                initial_relative_voltage_v,
                solve,
                cancelled,
            )?,
        })
    }

    /// Prescribed P in local (r,z), internal surfaces of revolution carrying
    /// free charge, and physical outward D.n traces use the same exact radial
    /// moments as the native linear field. The constitutive D(E) remains the
    /// original nonlinear law; fixed P contributes no inferred storage law.
    pub fn solve_with_sources(
        &self,
        sources: AxisymmetricElectrostaticSources<'_>,
        controls: &[ElectrodeControl],
        initial_relative_voltage_v: &[f64],
        solve: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricNonlinearElectricField<'m>, Error> {
        Ok(AxisymmetricNonlinearElectricField {
            problem: self.problem,
            field: self.core.solve_with_sources(
                PlanarElectrostaticSources {
                    free_charge_density_c_m3: sources.free_charge_density_c_m3,
                    impressed_polarization_c_m2: sources.impressed_polarization_c_m2,
                    charge_sheets: sources.charge_sheets,
                    displacement_traces: sources.displacement_traces,
                },
                controls,
                initial_relative_voltage_v,
                solve,
                cancelled,
            )?,
        })
    }
}

/// Accepted original nonlinear meridional field and physical conservation gates.
pub struct AxisymmetricNonlinearElectricField<'m> {
    problem: &'m AxisymmetricNonlinearElectrostatics,
    field: NonlinearElectricField<'m>,
}

impl<'m> AxisymmetricNonlinearElectricField<'m> {
    pub fn problem(&self) -> &'m AxisymmetricNonlinearElectrostatics {
        self.problem
    }
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.field.relative_voltage_v()
    }
    /// Physical reference attached separately to each original mesh vertex.
    pub fn reference_voltage_by_vertex_v(&self) -> &[f64] {
        self.field.reference_voltage_by_vertex_v()
    }
    pub fn reference_terminal_by_vertex(&self) -> &[usize] {
        self.field.reference_terminal_by_vertex()
    }
    pub fn reference_voltage_for_terminal_v(&self, terminal: usize) -> Option<f64> {
        self.field.reference_voltage_for_terminal_v(terminal)
    }
    /// Electrode potentials relative to their component's voltage reference.
    pub fn terminal_voltage_v(&self) -> &[f64] {
        self.field.terminal_voltage_v()
    }
    pub fn electrode_controls(&self) -> &[ElectrodeControl] {
        self.field.electrode_controls()
    }
    pub fn terminal_charge_c(&self) -> &[f64] {
        self.field.terminal_charge_c()
    }
    pub fn nodal_boundary_charge_c(&self) -> &[Option<f64>] {
        self.field.nodal_boundary_charge_c()
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        self.field.free_charge_density_c_m3()
    }
    /// Original cellwise `[E_r, E_z]` in V/m, with no azimuthal field.
    pub fn electric_field_v_m(&self) -> &[[f64; 2]] {
        self.field.electric_field_v_m()
    }
    /// Original cellwise `[D_r, D_z]` in C/m²; material jumps are retained.
    pub fn displacement_c_m2(&self) -> &[[f64; 2]] {
        self.field.displacement_c_m2()
    }
    pub fn constitutive_displacement_c_m2(&self) -> &[[f64; 2]] {
        self.field.constitutive_displacement_c_m2()
    }
    pub fn impressed_polarization_c_m2(&self) -> Option<&[[f64; 2]]> {
        self.field.impressed_polarization_c_m2()
    }
    pub fn free_charge_sheets(&self) -> &[ElectrostaticChargeSheet] {
        self.field.free_charge_sheets()
    }
    pub fn displacement_traces(&self) -> &[ElectrostaticDisplacementTrace] {
        self.field.displacement_traces()
    }
    pub fn sheet_charge_c(&self) -> f64 {
        self.field.sheet_charge_c()
    }
    pub fn outward_displacement_c(&self) -> f64 {
        self.field.outward_displacement_c()
    }
    pub fn polarization_pairing_j(&self) -> f64 {
        self.field.polarization_pairing_j()
    }
    pub fn cell_energy_j(&self) -> &[f64] {
        self.field.cell_energy_j()
    }
    pub fn cell_coenergy_j(&self) -> &[f64] {
        self.field.cell_coenergy_j()
    }
    /// Physical storage integral E.dD over the authored angular span [J].
    pub fn energy_j(&self) -> f64 {
        self.field.energy_j()
    }
    /// Original voltage potential integral D.dE; not additional storage [J].
    pub fn coenergy_j(&self) -> f64 {
        self.field.coenergy_j()
    }
    pub fn volume_charge_c(&self) -> f64 {
        self.field.volume_charge_c()
    }
    pub fn charge_balance_residual_c(&self) -> f64 {
        self.field.charge_balance_residual_c()
    }
    pub fn maximum_free_charge_residual_c(&self) -> f64 {
        self.field.maximum_free_charge_residual_c()
    }
    pub fn maximum_component_charge_residual_c(&self) -> f64 {
        self.field.maximum_component_charge_residual_c()
    }
    pub fn maximum_floating_charge_residual_c(&self) -> f64 {
        self.field.maximum_floating_charge_residual_c()
    }
    pub fn energy_identity_residual_j(&self) -> f64 {
        self.field.energy_identity_residual_j()
    }
    pub fn report(&self) -> &NonlinearElectrostaticReport {
        self.field.report()
    }

    /// Extract local dQ/dV at this operating point, with every floating total
    /// charge fixed. Reuses one prepared analytic radial tangent for all
    /// columns and retains the original reciprocity and charge checks.
    pub fn differential_capacitance(
        &self,
        controls: DifferentialCapacitanceControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricDifferentialCapacitance<'_, 'm>, Error> {
        Ok(AxisymmetricDifferentialCapacitance {
            field: self,
            response: self.field.differential_capacitance(controls, cancelled)?,
        })
    }

    /// First-order mechanical work [J] on a meridional mesh displacement.
    /// The angular span, voltage controls, floating charges and isotropic
    /// strain-independent material laws remain fixed. Free volume charge is
    /// refused: its moving radial P1 moments require an explicit transport law.
    ///
    /// This differentiates the actual discrete coenergy/source functional.
    /// Besides its planar gradient and area terms, it includes the variation
    /// of the radial Jacobian: W*_cell * mean(delta r)/mean(r). Omitting this
    /// term gives incorrect radial force even in a uniform axial field.
    /// The result equals minus the storage derivative under fixed charge.
    /// Axis vertices must remain on the axis. This is a virtual derivative,
    /// not a finite trajectory, nodal solid load or continuum force certificate.
    pub fn virtual_work(
        &self,
        displacement_m: &[[f64; 2]],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        poll(cancelled)?;
        let mesh = self.problem.mesh();
        if self.field.additional_sources.has_geometry_sources() {
            return Err(invalid(
                "axisymmetric shape work with impressed polarization or surface sources requires their transport law",
            ));
        }
        if displacement_m.len() != mesh.vertices().len()
            || displacement_m.iter().flatten().any(|v| !v.is_finite())
            || self
                .field
                .free_charge_density_c_m3
                .iter()
                .any(|&rho| rho != 0.0)
            || mesh
                .vertices()
                .iter()
                .zip(displacement_m)
                .any(|(x, u)| x[0] == 0.0 && u[0] != 0.0)
        {
            return Err(invalid(
                "axisymmetric virtual work needs finite meridional motion, a fixed axis and no transported volume charge",
            ));
        }
        let (mut work, mut correction) = (0.0, 0.0);
        for (cell, element) in self.problem.core.elements.iter().enumerate() {
            poll(cancelled)?;
            let forces = self.field.shape_cell(cell)?;
            let u = element.vertices.map(|v| displacement_m[v]);
            for local in 0..2 {
                for axis in 0..2 {
                    add(
                        &mut work,
                        &mut correction,
                        forces[local][axis] * finite(u[local + 1][axis] - u[0][axis])?,
                    )?;
                }
            }
            // Scaling the three positive radii before summation avoids
            // overflowing their sum. The same rule bounds the displacement
            // mean before its signed compensated reduction.
            let mut radius = 0.0;
            let (mut radial_motion, mut motion_correction) = (0.0, 0.0);
            for (&vertex, motion) in element.vertices.iter().zip(u) {
                radius += mesh.vertices()[vertex][0] / 3.0;
                add(&mut radial_motion, &mut motion_correction, motion[0] / 3.0)?;
            }
            let radial_ratio = finite(finite(radial_motion + motion_correction)? / radius)?;
            add(
                &mut work,
                &mut correction,
                self.field.evaluation.cell_coenergy[cell] * radial_ratio,
            )?;
        }
        poll(cancelled)?;
        Ok(finite(work + correction)?)
    }
}

/// Source-bound local differential terminal response for a revolved field.
/// Its field accessor preserves the axisymmetric geometry boundary.
pub struct AxisymmetricDifferentialCapacitance<'a, 'm> {
    field: &'a AxisymmetricNonlinearElectricField<'m>,
    response: DifferentialCapacitance<'a, 'm>,
}

impl<'a, 'm> AxisymmetricDifferentialCapacitance<'a, 'm> {
    pub fn field(&self) -> &'a AxisymmetricNonlinearElectricField<'m> {
        self.field
    }
    pub fn terminals(&self) -> &[usize] {
        self.response.terminals()
    }
    pub fn matrix_f(&self) -> &[f64] {
        self.response.matrix_f()
    }
    pub fn maximum_residual_f(&self) -> f64 {
        self.response.maximum_residual_f()
    }
    pub fn operator_applications(&self) -> usize {
        self.response.operator_applications()
    }
}
