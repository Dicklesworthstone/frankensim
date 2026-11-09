//! Original measured-source liquid expansion on the native Stokes constraint.

use super::{CaloricDilationError, StokesAssembly, StokesSolution, solve_dilation};
use crate::ns::boundary::StokesBoundary;
use crate::ns::{FluxParams, FluxSystem};
use fs_material::thermal::transport::{
    ThermalTransportError,
    isobaric::expansion::liquid::{LiquidIsobaricDilation, LiquidIsobaricExpansion},
};

/// One original source isobar and its material enthalpy rate per flow cell.
#[derive(Clone, Copy)]
pub struct IsobaricDilationCell<'a> {
    pub caloric: &'a LiquidIsobaricExpansion<'a>,
    pub specific_enthalpy_j_kg: f64,
    pub material_specific_enthalpy_rate_w_kg: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricDilationError {
    Flow(CaloricDilationError),
    Material {
        cell: usize,
        source: ThermalTransportError,
    },
}
impl core::fmt::Display for IsobaricDilationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original-source Stokes dilation: {self:?}")
    }
}
impl std::error::Error for IsobaricDilationError {}
impl From<CaloricDilationError> for IsobaricDilationError {
    fn from(error: CaloricDilationError) -> Self {
        Self::Flow(error)
    }
}

/// Original native matrix, with source-derived integrated continuity loads.
pub struct IsobaricDilationAssembly<'a, 'm> {
    native: StokesAssembly<'a, 'm>,
    source: Vec<LiquidIsobaricDilation>,
}

/// Native velocity/pressure and original source tangents used in every cell.
pub struct IsobaricDilationSolution {
    pub flow: StokesSolution,
    pub source: Vec<LiquidIsobaricDilation>,
    pub maximum_dilation_defect_per_s: f64,
}

impl<'a, 'm> IsobaricDilationAssembly<'a, 'm> {
    #[must_use]
    pub const fn native(&self) -> &StokesAssembly<'a, 'm> {
        &self.native
    }
    #[must_use]
    pub fn source(&self) -> &[LiquidIsobaricDilation] {
        &self.source
    }

    /// Solve and independently read every actual BDM cell flux, using the same
    /// physical allowance and component/gauge treatment as caloric charts.
    pub fn solve(
        &self,
        maximum_dilation_defect_per_s: f64,
    ) -> Result<IsobaricDilationSolution, IsobaricDilationError> {
        let (flow, maximum_defect) = solve_dilation(
            &self.native,
            self.source.iter().map(|source| source.dilation_per_s),
            maximum_dilation_defect_per_s,
        )?;
        Ok(IsobaricDilationSolution {
            flow,
            source: self.source.clone(),
            maximum_dilation_defect_per_s: maximum_defect,
        })
    }
}

impl<'m> FluxSystem<'m> {
    /// Instantaneous liquid flow from the ORIGINAL Cp and engineering-length
    /// sources. The admitted material path supplies rho=rho_ref/J and exact
    /// directional dv/dh, including left-sided cooling at genuine source knots.
    /// No replacement density/CTE fit or resampled phase chart is constructed.
    ///
    /// The symmetric-gradient native operator, stress/effective-pressure
    /// convention, traction outlets, sealed-component compatibility and fixed
    /// reference-density normalization are exactly those documented by
    /// `assemble_stokes_caloric_dilation`. Thermodynamic source pressure stays
    /// fixed; native mechanical pressure is not a new compression EOS. The
    /// entry requires one identical source pressure per native connected
    /// component; disconnected components may have different isobars. No
    /// material pressure-derivative closure is inferred from mechanical p.
    /// This snapshot advances neither mass, enthalpy, geometry nor material transport.
    pub fn assemble_stokes_isobaric_dilation<'a>(
        &'a self,
        params: FluxParams,
        viscosity: &dyn Fn(usize, [f64; 2]) -> f64,
        force: &dyn Fn([f64; 2]) -> [f64; 2],
        boundaries: &[StokesBoundary<'_>],
        cells: &[IsobaricDilationCell<'_>],
    ) -> Result<IsobaricDilationAssembly<'a, 'm>, IsobaricDilationError> {
        if cells.len() != self.mesh.tris.len() {
            return Err(CaloricDilationError::Invalid(
                "one original isobaric source and material heat rate per triangle required",
            )
            .into());
        }
        let mut source = Vec::with_capacity(cells.len());
        let mut outflow = Vec::with_capacity(cells.len());
        let mut pressures = vec![None; self.pressure_roots.len()];
        for (cell, input) in cells.iter().enumerate() {
            let pressure = input.caloric.caloric().caloric().pressure_pa();
            let component = self.pressure_component[cell];
            if pressures[component].is_some_and(|original| original != pressure) {
                return Err(CaloricDilationError::Invalid(
                    "connected liquid cells require one original maintained source pressure",
                )
                .into());
            }
            pressures[component] = Some(pressure);
            let response = input
                .caloric
                .dilation_at_specific_enthalpy(
                    input.specific_enthalpy_j_kg,
                    input.material_specific_enthalpy_rate_w_kg,
                )
                .map_err(|source| IsobaricDilationError::Material { cell, source })?;
            let flux = self.mesh.areas[cell] * response.dilation_per_s;
            if !flux.is_finite() || (response.dilation_per_s != 0.0 && flux == 0.0) {
                return Err(
                    CaloricDilationError::Invalid("unrepresentable thermal outflow").into(),
                );
            }
            source.push(response);
            outflow.push(flux);
        }
        let native =
            self.assemble_stokes_dilation_native(params, viscosity, force, boundaries, &outflow)?;
        Ok(IsobaricDilationAssembly { native, source })
    }
}
