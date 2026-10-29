//! Prescribed gas changes in the continuing valve and original DG owners.
//!
//! Native gas transfer preserves compression and axial momentum, so acoustic
//! pressure and volume flow change. The original DG owner measures the storage
//! change as signed external parameter work. The quasisteady valve carries no
//! gas inertia: updating its density, source viscosity and actual mean impedance
//! changes no mechanical/contact/material coordinate or stored energy.

use super::{ApertureDg, ApertureDgLoad, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::two_sided::LoadedAperture;
use crate::dg_port::DgGasRefresh;
use fs_duct::dg::DgGasMaterialUpdate;
use fs_material::gas_field::AxialGasField;

fn prepare_single<A>(
    model: &ApertureDg<A>,
    field: &AxialGasField,
    refresh: DgGasRefresh,
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ApertureDgLoad, DgGasMaterialUpdate), AcousticRealizeError> {
    if model.upstream.is_some() {
        return Err(invalid(
            "live two-sided DG gas refresh requires both explicit gas fields in one transaction",
        ));
    }
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    let mut load = model.load.clone();
    let report = load.rebind_gas_field_controlled(
        field,
        refresh,
        maximum_absolute_parameter_work_j,
        cancelled,
    )?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    Ok((load, report))
}

impl ApertureDg {
    /// Refresh the running bore and moving valve from one complete gas field.
    ///
    /// The original native transfer retains compression, momentum, wall and
    /// diffusion histories and records its actual storage change as parameter
    /// work. The source uses the SAME new inlet density and exact mean-port
    /// impedance. A source-bound viscous slit follows the new gas viscosity;
    /// an explicitly supplied viscosity remains the caller's explicit override.
    /// Mechanical/contact/solid-relaxation histories and sample clocks survive.
    ///
    /// Rational exterior and independent branch gases remain their own original
    /// sources. Their histories are retained, and the complete mean-port relation
    /// is recomputed. This prescribed intervention does not solve gas heating,
    /// species transport, static mean-pressure traction or radiated timbre.
    /// Every fallible check precedes publication of either owner.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        refresh: DgGasRefresh,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            field,
            refresh,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        // This update validates before touching the original source. There is
        // no fallible operation between source and complete load publication.
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.load = load;
        Ok(report)
    }
}

impl ApertureDg<LoadedAperture> {
    /// The same atomic live gas transfer for original single, independent-pair
    /// or material-modal mechanics. Full source viscosity, contact, signed mode
    /// projections and solid constitutive memories stay in their original owner.
    ///
    /// With an upstream tract, use `rebind_series_gas_fields_controlled` so both
    /// complete inlet gases and the summed load are admitted together.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        refresh: DgGasRefresh,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            field,
            refresh,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.load = load;
        Ok(report)
    }

    /// Refresh both continuing DG fields and their common valve atomically.
    ///
    /// Inputs and returned reports are ordered UPSTREAM, DOWNSTREAM. Each field
    /// covers its own original geometry; remote composition and loss models may
    /// differ. Complete inlet gases must match, since the existing small-signal
    /// volume-flow junction has no density-changing mass-transport model. The
    /// valve uses the SUM of the two new exact mean impedances.
    ///
    /// The finite nonnegative budget bounds the SUM OF ABSOLUTE native parameter
    /// work on both sides: opposite interventions cannot cancel an exceeded
    /// allowance. The signed reports enter their own original ledgers once.
    /// All accepted acoustic, branch, exterior, wall and mechanical histories
    /// survive. Failure or cancellation publishes neither gas field nor source.
    pub fn rebind_series_gas_fields_controlled(
        &mut self,
        fields: [&AxialGasField; 2],
        refresh: [DgGasRefresh; 2],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[DgGasMaterialUpdate; 2], AcousticRealizeError> {
        let Some(upstream) = &self.upstream else {
            return Err(invalid(
                "live series gas refresh requires an original upstream DG tract",
            ));
        };
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let mut up = upstream.clone();
        let mut down = self.load.clone();
        let up_report = up.rebind_gas_field_controlled(
            fields[0],
            refresh[0],
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let down_report = down.rebind_gas_field_controlled(
            fields[1],
            refresh[1],
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let gas = *down.duct().inlet_gas();
        let impedance = up.impedance_pa_s_m3() + down.impedance_pa_s_m3();
        let absolute_work = up_report.parameter_work_j.abs() + down_report.parameter_work_j.abs();
        if up.duct().inlet_gas() != &gas
            || !impedance.is_finite()
            || impedance <= 0.0
            || !absolute_work.is_finite()
            || absolute_work > maximum_absolute_parameter_work_j
        {
            return Err(invalid(
                "live DG tract refresh requires one complete inlet gas, finite summed impedance and admitted total absolute parameter work",
            ));
        }
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        self.aperture.rebind_exact_gas_load(&gas, impedance)?;
        self.upstream = Some(up);
        self.load = down;
        Ok([up_report, down_report])
    }
}
