//! Resolved circular-bore fields with actual finite/layered thermal walls.
//!
//! The disk has one physical boundary and a genuine regular axis at r=0.
//! It uses the SAME exact r-weighted P1 mass/stiffness assembly, complementary
//! diffusion variables, native tridiagonal solve and gas/molecular energy
//! equation as the original annular owner. The axis remains a free coordinate:
//! no fictitious inner no-slip/isothermal wall or small-hole limit is introduced.
//!
//! The declared uniform radial mesh is the spatial accuracy axis. Source-backed
//! finite-thickness and layered solids retain their independently prepared
//! radial meshes, actual material coefficients, contact losses and bath.
//! A material wall must touch the circular gas at its inner face, at the same
//! radius, frequency and background temperature. The gas's original dynamic
//! Cp enters once. Wall coefficients and acoustic loss are never fitted.
//!
//! This is linear, axisymmetric, fully developed, no-slip, plane-pressure
//! continuum acoustics. The mechanical wall remains rigid. End fields, axial
//! heat diffusion, mean flow, rarefied slip, nonuniform background temperature
//! and thermoelastic motion require other owners. Residual and phase admission
//! do not replace radial refinement or establish empirical material accuracy.

use super::annular::{
    AnnularSection,
    thermal_wall::conduction::{
        PreparedConjugateAnnularWalls, RadialThermalWallError, RadialThermalWallField,
        layered::RadialThermalContactObservation,
    },
};
use crate::{DuctError, SegmentWave, compose::TwoPort};
use fs_material::{gas::GasState, relaxation::RelaxationState};
use fs_math::c64::C64;

/// Original explicit plane-pressure, radial resolution, residual and byte bounds.
/// The existing conservative transverse bound uses the actual outer radius.
pub use super::annular::AnnularAdmission as CircularAdmission;
/// Explicit isothermal, insulated or finite-capacity bath boundary.
pub use super::annular::thermal_wall::AnnularThermalWall as CircularThermalBoundary;
/// Actual gas velocity/temperature and gradients, with native peak phasors.
pub use super::annular::thermal_wall::AnnularThermalWallField as CircularThermalField;
/// Original explicit boundary, homogeneous solid, or layered material wall.
pub use super::annular::thermal_wall::conduction::AnnularConjugateWall as CircularThermalWall;
/// One physical wall's actual face fluxes, losses and reversible thermal storage.
pub use super::annular::thermal_wall::conduction::AnnularConjugateWallObservation as CircularThermalWallObservation;

/// Prepared disk geometry with uniform radial P1 cells, retaining the axis.
#[derive(Clone, Debug)]
pub struct CircularSection {
    radial: AnnularSection,
}

impl CircularSection {
    pub fn new(radius_m: f64, elements: usize) -> Result<Self, DuctError> {
        Self::new_controlled(radius_m, elements, || false)
    }

    /// Positive radius and 2..=4096 original radial elements. Geometry work
    /// polls before allocation and every 64 cells; no partial section escapes.
    pub fn new_controlled(
        radius_m: f64,
        elements: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DuctError> {
        Ok(Self {
            radial: AnnularSection::circular_controlled(radius_m, elements, cancelled)?,
        })
    }

    pub fn radius_m(&self) -> f64 {
        self.radial.outer_radius_m()
    }
    pub fn area_m2(&self) -> f64 {
        self.radial.area_m2()
    }
    pub fn elements(&self) -> usize {
        self.radial.elements()
    }
    /// Physical radius of an original radial node, including zero at the axis.
    pub fn node_radius_m(&self, node: usize) -> Option<f64> {
        self.radial.radius_m(node)
    }

    pub fn prepare_thermal_wall_fields<'g, 'm, 'w>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        wall: CircularThermalWall<'w>,
        admission: CircularAdmission,
    ) -> Result<PreparedCircularThermalWalls<'g, 'm, 'w, GasState>, RadialThermalWallError> {
        self.prepare_thermal_wall_fields_controlled(gas, omega, wall, admission, || false)
    }

    /// The shared O(N) gas solve requires 320*(elements+1)+4096 payload bytes.
    /// Separately prepared solid/material owners retain their own budgets.
    /// The regular axis has zero physical boundary measure and no heat loss.
    pub fn prepare_thermal_wall_fields_controlled<'g, 'm, 'w>(
        &'g self,
        gas: &'m GasState,
        omega: f64,
        wall: CircularThermalWall<'w>,
        admission: CircularAdmission,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCircularThermalWalls<'g, 'm, 'w, GasState>, RadialThermalWallError> {
        let fields = self.radial.prepare_conjugate_wall_fields_controlled(
            gas,
            omega,
            [
                CircularThermalWall::Boundary(CircularThermalBoundary::Adiabatic),
                wall,
            ],
            admission,
            cancelled,
        )?;
        Ok(PreparedCircularThermalWalls {
            section: self,
            fields,
        })
    }

    pub fn prepare_molecular_thermal_wall_fields<'g, 'm, 'w>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        wall: CircularThermalWall<'w>,
        admission: CircularAdmission,
    ) -> Result<PreparedCircularThermalWalls<'g, 'm, 'w, RelaxationState>, RadialThermalWallError>
    {
        self.prepare_molecular_thermal_wall_fields_controlled(
            material,
            omega,
            wall,
            admission,
            || false,
        )
    }

    /// The same regular-axis shear solve and thermal equation, with the
    /// original molecular Cp and existing solid response coupled exactly once.
    pub fn prepare_molecular_thermal_wall_fields_controlled<'g, 'm, 'w>(
        &'g self,
        material: &'m RelaxationState,
        omega: f64,
        wall: CircularThermalWall<'w>,
        admission: CircularAdmission,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCircularThermalWalls<'g, 'm, 'w, RelaxationState>, RadialThermalWallError>
    {
        let fields = self
            .radial
            .prepare_molecular_conjugate_wall_fields_controlled(
                material,
                omega,
                [
                    CircularThermalWall::Boundary(CircularThermalBoundary::Adiabatic),
                    wall,
                ],
                admission,
                cancelled,
            )?;
        Ok(PreparedCircularThermalWalls {
            section: self,
            fields,
        })
    }
}

/// Independent positive loss channels per axial metre. The axis contributes
/// none; thermal-contact loss belongs only to explicit interlayer resistance.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CircularThermalWallDissipation {
    pub gas_viscous_w_m: f64,
    pub gas_conduction_w_m: f64,
    pub gas_molecular_w_m: f64,
    pub solid_conduction_w_m: f64,
    pub thermal_contact_w_m: f64,
    pub bath_exchange_w_m: f64,
    pub total_w_m: f64,
}

/// Original gas and material wall, bound to one actual circular radial mesh.
/// Local, wall and integrated queries allocate nothing and are O(1).
#[derive(Debug)]
pub struct PreparedCircularThermalWalls<'g, 'm, 'w, M: ?Sized> {
    section: &'g CircularSection,
    fields: PreparedConjugateAnnularWalls<'g, 'm, 'w, M>,
}

impl<'g, 'm, 'w, M: ?Sized> PreparedCircularThermalWalls<'g, 'm, 'w, M> {
    pub fn source(&self) -> &'m M {
        self.fields.source()
    }
    pub const fn section(&self) -> &'g CircularSection {
        self.section
    }
    pub fn wall(&self) -> &CircularThermalWall<'w> {
        &self.fields.walls()[1]
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.fields.angular_frequency_rad_s()
    }
    pub fn wave(&self) -> &SegmentWave {
        self.fields.wave()
    }
    /// Area-mean normalized thermal response, in +i omega solve convention.
    pub fn mean_temperature_response(&self) -> C64 {
        self.fields.mean_temperature_response()
    }
    pub fn relative_backward_residual(&self) -> f64 {
        self.fields.relative_backward_residual()
    }
    pub fn relative_thermal_work_defect(&self) -> f64 {
        self.fields.relative_thermal_work_defect()
    }

    pub fn two_port(&self, length_m: f64, reference: f64) -> Result<TwoPort, DuctError> {
        self.fields.two_port(length_m, reference)
    }

    /// The original P1 radial field. Axis values are unconstrained gas values,
    /// while adjacent-cell derivatives, including the first cell, converge
    /// with radial refinement; no exact pointwise derivative is fabricated.
    pub fn field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<CircularThermalField, DuctError> {
        self.fields
            .field_in_cell(cell, fraction, pressure_pa, volume_flow_m3_s)
    }

    pub fn wall_observation(
        &self,
        pressure_pa: C64,
    ) -> Result<CircularThermalWallObservation, DuctError> {
        self.fields.wall_observation(1, pressure_pa)
    }

    pub fn solid_field_in_cell(
        &self,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<RadialThermalWallField, DuctError> {
        self.fields
            .solid_field_in_cell(1, cell, fraction, pressure_pa)
    }

    pub fn layered_solid_field_in_cell(
        &self,
        layer: usize,
        cell: usize,
        fraction: f64,
        pressure_pa: C64,
    ) -> Result<RadialThermalWallField, DuctError> {
        self.fields
            .layered_solid_field_in_cell(1, layer, cell, fraction, pressure_pa)
    }

    pub fn thermal_contact_observation(
        &self,
        contact: usize,
        pressure_pa: C64,
    ) -> Result<RadialThermalContactObservation, DuctError> {
        self.fields
            .thermal_contact_observation(1, contact, pressure_pa)
    }

    pub fn dissipation_per_axial_length(
        &self,
        pressure_pa: C64,
        volume_flow_m3_s: C64,
    ) -> Result<CircularThermalWallDissipation, DuctError> {
        let original = self
            .fields
            .dissipation_per_axial_length(pressure_pa, volume_flow_m3_s)?;
        Ok(CircularThermalWallDissipation {
            gas_viscous_w_m: original.gas_viscous_w_m,
            gas_conduction_w_m: original.gas_conduction_w_m,
            gas_molecular_w_m: original.gas_molecular_w_m,
            solid_conduction_w_m: original.solid_conduction_w_m[1],
            thermal_contact_w_m: original.thermal_contact_w_m[1],
            bath_exchange_w_m: original.bath_exchange_w_m[1],
            total_w_m: original.total_w_m,
        })
    }
}
