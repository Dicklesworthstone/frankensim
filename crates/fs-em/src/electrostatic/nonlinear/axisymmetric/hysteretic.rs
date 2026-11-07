//! Retained dielectric polarization on the original revolved P1 complex.
//!
//! The immutable isotropic source law acts in local cylindrical coordinates.
//! Native three-vector histories are ordered `[r, z, theta]`; their third
//! component stays zero. A radial polarization therefore follows the radial
//! basis around the authored angular span, rather than a fixed Cartesian axis.
//! Geometry, sources and that interpretation remain fixed throughout a path.
//!
//! The shared Gauss and vector-play owners retain the exact measure
//! `angle * r * dr * dz`, physical switching heat, endpoint source work and a
//! separate numerical loss. Changing volume charge uses the same radial P1
//! moments as the field equation. Angular spans are symmetry sectors, with no
//! independently driven azimuthal faces. Rates, leakage, material-coefficient
//! changes, electrocaloric response and moving geometry are not inferred.

pub use super::super::hysteretic::{
    DielectricSwitchingCell, HystereticElectrostaticStep,
    path::{
        HystereticPathControls, HystereticPathError, HystereticPathReport, HystereticPathSummary,
    },
};
use super::super::{
    DielectricMaterial, ElectrodeControl, Error, NonlinearElectrostaticControls,
    hysteretic::PlanarHystereticElectrostatics, invalid, poll,
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::QueryPoint;
use fs_material::field_storage::dielectric::hysteresis::vector::{
    ResolvedVectorDielectricHysteresis, VectorDielectricHysteresisControls,
    VectorDielectricHysteresisState,
};
use fs_rep_mesh::{Metric2, TriComplex2};

/// Accepted meridional field and polarization memory over the authored angular
/// span. The shared implementation remains private: the public field cannot be
/// passed to planar force or planar caloric coupling owners.
#[derive(Clone)]
pub struct AxisymmetricHystereticElectrostatics {
    core: PlanarHystereticElectrostatics,
}

impl AxisymmetricHystereticElectrostatics {
    /// Start at virgin zero field, free charge and polarization. Every original
    /// isotropic fit is interpreted in the local `(r,z)` material directions.
    /// Terminal traces must have positive revolved surface measure; an axis
    /// line alone cannot become a voltage or total-charge electrode.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedVectorDielectricHysteresis>,
        terminals: &[&[usize]],
        local: VectorDielectricHysteresisControls,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(cancelled)?;
        if !matches!(mesh.metric(), Metric2::Axisymmetric { .. }) {
            return Err(invalid(
                "cylindrical retained polarization requires an explicit radial metric",
            ));
        }
        Ok(Self {
            core: PlanarHystereticElectrostatics::assemble_metric(
                mesh, materials, terminals, local, limits, cancelled,
            )?,
        })
    }

    /// Original metre-valued `[radius, axial_coordinate]` geometry.
    pub fn mesh(&self) -> &TriComplex2 {
        self.core.problem().mesh()
    }
    /// Original isotropic source laws and their accepted local histories.
    pub fn materials(&self) -> &[DielectricMaterial] {
        self.core.problem().materials()
    }
    pub fn terminals(&self) -> &[Vec<usize>] {
        self.core.terminals()
    }
    /// Potentials relative to each vertex's electrical-component reference.
    pub fn relative_voltage_v(&self) -> &[f64] {
        self.core.relative_voltage_v()
    }
    pub fn reference_terminal_by_vertex(&self) -> &[usize] {
        self.core.reference_terminal_by_vertex()
    }
    pub fn reference_voltage_by_vertex_v(&self) -> &[f64] {
        self.core.reference_voltage_by_vertex_v()
    }
    pub fn reference_terminal_for_terminal(&self, terminal: usize) -> Option<usize> {
        self.core.reference_terminal_for_terminal(terminal)
    }
    pub fn reference_voltage_for_terminal_v(&self, terminal: usize) -> Option<f64> {
        self.core.reference_voltage_for_terminal_v(terminal)
    }
    /// Each electrode's voltage relative to its own component's reference.
    pub fn terminal_voltage_v(&self) -> &[f64] {
        self.core.terminal_voltage_v()
    }
    /// Total conductor free charge over the original angular span [C].
    pub fn terminal_charge_c(&self) -> &[f64] {
        self.core.terminal_charge_c()
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.core.stored_energy_j()
    }
    pub fn free_charge_density_c_m3(&self) -> &[f64] {
        self.core.free_charge_density_c_m3()
    }
    /// Accepted E, D and every branch polarization, in local `[r,z,theta]`
    /// order. The azimuthal component is zero, including on the symmetry axis.
    pub fn cell_history(&self, triangle: usize) -> Option<&VectorDielectricHysteresisState> {
        self.core.cell_history(triangle)
    }

    /// Recheck the same constant fits at new material queries without changing
    /// geometry, fields, charges, storage or polarization. Changed coefficients
    /// refuse through the original source adapter; this is not thermal evolution.
    pub fn at_source_queries(
        &self,
        queries: &[QueryPoint],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Ok(Self {
            core: self.core.at_source_queries(queries, cancelled)?,
        })
    }

    /// Advance one bounded quasistatic endpoint with original voltage or TOTAL
    /// free-charge controls. Sources are per physical revolved cell volume.
    /// One voltage reference is required per electrical component. Refusal
    /// leaves every accepted cell and independent gauge unchanged.
    pub fn solve_candidate(
        &self,
        free_charge_density_c_m3: &[f64],
        boundary: &[ElectrodeControl],
        controls: NonlinearElectrostaticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricHystereticElectrostaticCandidate, Error> {
        let (core, step) = self
            .core
            .solve_candidate(free_charge_density_c_m3, boundary, controls, cancelled)?
            .into_parts();
        Ok(AxisymmetricHystereticElectrostaticCandidate {
            field: Self { core },
            step,
        })
    }

    /// Resolve one affine loading segment using the original coarse/two-half
    /// error estimates and shared solve/operator budgets. Electrical control
    /// modes stay fixed within a segment; separate segments describe reversals
    /// and mode changes. Only accepted substeps become physical history.
    pub fn solve_affine_path_candidate(
        &self,
        target_free_charge_density_c_m3: &[f64],
        target_controls: &[ElectrodeControl],
        controls: HystereticPathControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricHystereticPathCandidate, HystereticPathError> {
        let (core, report) = self
            .core
            .solve_affine_path_candidate(
                target_free_charge_density_c_m3,
                target_controls,
                controls,
                cancelled,
            )?
            .into_parts();
        Ok(AxisymmetricHystereticPathCandidate {
            field: Self { core },
            report,
        })
    }
}

/// Complete cylindrical endpoint and its independent physical work/heat account.
pub struct AxisymmetricHystereticElectrostaticCandidate {
    field: AxisymmetricHystereticElectrostatics,
    step: HystereticElectrostaticStep,
}

impl AxisymmetricHystereticElectrostaticCandidate {
    pub fn field(&self) -> &AxisymmetricHystereticElectrostatics {
        &self.field
    }
    pub fn step(&self) -> &HystereticElectrostaticStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        AxisymmetricHystereticElectrostatics,
        HystereticElectrostaticStep,
    ) {
        (self.field, self.step)
    }
}

/// Accepted cylindrical path; the report sums physical switching heat and
/// numerical loss separately over the actual revolved cells.
pub struct AxisymmetricHystereticPathCandidate {
    field: AxisymmetricHystereticElectrostatics,
    report: HystereticPathReport,
}

impl AxisymmetricHystereticPathCandidate {
    pub fn field(&self) -> &AxisymmetricHystereticElectrostatics {
        &self.field
    }
    pub fn report(&self) -> &HystereticPathReport {
        &self.report
    }
    pub fn into_parts(self) -> (AxisymmetricHystereticElectrostatics, HystereticPathReport) {
        (self.field, self.report)
    }
    /// Retain the cumulative native account while keeping cylindrical ownership.
    /// This does not coerce the field into a planar heat or force adapter.
    pub fn into_electrostatic_candidate(
        self,
    ) -> (
        AxisymmetricHystereticElectrostaticCandidate,
        HystereticPathSummary,
    ) {
        (
            AxisymmetricHystereticElectrostaticCandidate {
                field: self.field,
                step: self.report.account,
            },
            self.report.summary,
        )
    }
}
