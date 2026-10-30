//! Complete shared heat transport on the original planar induction history.
use super::*;
use crate::circuit_enthalpy::{
    CircuitEnthalpyError,
    contacts::heat::{HeatNetworkError, HeatNetworkSample, PreparedHeatNetwork, closed_balance},
    eddy::nonlinear::thermal::NonlinearEddyCaloricCandidate,
    magnetic::{MagneticRegion, MagneticThermalError, caloric::MagneticConstitutiveThermalLaw},
};
use fs_em::planar::{
    eddy::nonlinear::{NonlinearEddyControls, PlanarNonlinearEddy},
    nonlinear::IsotropicMagneticMaterial,
};
use fs_exec::Cx;

pub(in crate::circuit_enthalpy) fn heat_error(error: HeatNetworkError) -> Error {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MagneticThermalError::Heat(other)),
    }
}

pub(in crate::circuit_enthalpy) fn close_heat_balance<D>(
    report: &mut EddyCaloricReport<HeatNetworkSample, D>, dt: f64, tolerance: f64,
) -> Result<(), Error> {
    let sample = &report.coupled.source.material.network;
    closed_balance(dt, report.coupled.total_energy_residual_j, sample, tolerance)
        .map_err(heat_error)?;
    closed_balance(dt, report.total_energy_residual_j, sample, tolerance)
        .map_err(heat_error)?;
    let internal_net = finite(dt * sample.radiation_net_power_w())?;
    report.external_heat_j = finite(report.external_heat_j - internal_net)?;
    report.total_energy_residual_j = finite(report.total_energy_residual_j + internal_net)?;
    Ok(())
}

impl PreparedHeatNetwork<'_, '_> {
    /// Continue the SAME planar field and finite heat history with sourced
    /// thermal interfaces, integrated layer conduction and enclosure radiation.
    /// Their coefficients use the same trial temperatures as sigma(T) and the
    /// magnetic entropy law. Native Joule heat is deposited exactly once.
    ///
    /// The retained material assembly owns the caloric cells and all links;
    /// supplied regions bind those original volumes to the field triangles.
    /// Graph electrical ports refuse because this field step does not own their
    /// currents/work. Thermal boundaries must not depend on unowned E/J/B/H.
    /// Geometry, contact areas and view factors remain fixed. The existing
    /// first-order temporal coupling and source validity limits still apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_eddy_thermal_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarEddy<'m>,
        impressed_jz_a_m2: &[f64],
        regions: &[MagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[MagneticThermalLaw],
        controls: CircuitEnthalpySolve,
        field_controls: EddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
    ) -> Result<EddyCaloricCandidate<'t, 'm, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        self.admit_field_independent_boundaries(cx, &[]).map_err(heat_error)?;
        let layout = MagneticThermalNetwork {
            cells: bulk.cells, links: bulk.links, magnetic_regions: regions,
            volume_tolerance_m3,
        };
        let mut next = layout.step_eddy_thermal_resolved_candidate(
            thermal, target_s, field, impressed_jz_a_m2, bindings, controls,
            field_controls, assembly, initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid("planar induction cannot discard graph electrical work"));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        close_heat_balance(&mut next.report, target_s - thermal.time_s(),
            controls.total_energy_tolerance_j)?;
        check(&cancel)?;
        Ok(next)
    }

    /// Continue a mixed linear/nonlinear magnetic assembly with the SAME
    /// complete finite heat network. Each region declares its sourced thermal
    /// constitutive law; the existing nonlinear field owner preserves its
    /// retained potential/current history and reports numerical loss separately
    /// from Joule heat. A separable reference curve is an explicit material
    /// model, never an inferred arbitrary H(B,T) surface or a fitted core loss.
    #[allow(clippy::too_many_arguments)]
    pub fn step_nonlinear_eddy_thermal_candidate<'t, 'm>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        field: &PlanarNonlinearEddy<'m>,
        impressed_jz_a_m2: &[f64],
        regions: &[MagneticRegion],
        volume_tolerance_m3: f64,
        bindings: &[MagneticThermalLaw<MagneticConstitutiveThermalLaw>],
        controls: CircuitEnthalpySolve,
        field_controls: NonlinearEddyControls,
        assembly: WeightedAssemblyLimits,
        initial_energy_increments_j: &[f64],
    ) -> Result<NonlinearEddyCaloricCandidate<'t, 'm, HeatNetworkSample>, Error> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        check(&cancel)?;
        let bulk = self.thermal_material_network().map_err(heat_error)?;
        let mut magnetic_axes = Vec::new();
        for region in regions {
            check(&cancel)?;
            if let MagneticRegion::Cell { triangle, cell } = *region {
                if let Some(IsotropicMagneticMaterial::Reversible(law)) =
                    field.problem().materials().get(triangle) {
                    magnetic_axes.push((cell, law.source_path().axis.as_str()));
                }
            }
        }
        self.admit_field_independent_boundaries(cx, &magnetic_axes).map_err(heat_error)?;
        let layout = MagneticThermalNetwork {
            cells: bulk.cells, links: bulk.links, magnetic_regions: regions,
            volume_tolerance_m3,
        };
        let mut next = layout.step_nonlinear_eddy_thermal_resolved_candidate(
            thermal, target_s, field, impressed_jz_a_m2, bindings, controls,
            field_controls, assembly, initial_energy_increments_j,
            |trial| {
                let coefficients = self.resolve(cx, trial).map_err(heat_error)?;
                if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
                    return Err(invalid("nonlinear induction cannot discard graph electrical work"));
                }
                Ok(coefficients.bulk)
            },
            |sample| &sample.material,
            &cancel,
        )?;
        close_heat_balance(&mut next.report, target_s - thermal.time_s(),
            controls.total_energy_tolerance_j)?;
        check(&cancel)?;
        Ok(next)
    }
}
