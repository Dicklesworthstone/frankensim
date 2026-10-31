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
use crate::dg_port::{DgBranchGasInput, DgGasRefresh, DgPortGasUpdate, DgRegionalGasUpdate};
use fs_material::gas::GasState;
use fs_material::gas_field::AxialGasField;

fn prepare_single<A, R, C: FnMut() -> bool>(
    model: &ApertureDg<A>,
    update: impl FnOnce(&mut ApertureDgLoad, &mut C) -> Result<R, AcousticRealizeError>,
    cancelled: &mut C,
) -> Result<(ApertureDgLoad, R), AcousticRealizeError> {
    if model.supply_channel.is_some() {
        return Err(invalid("live DG gas refresh must include the fixed supply passage's momentum and parameter work"));
    }
    if model.upstream.is_some() {
        return Err(invalid(
            "live two-sided DG gas refresh requires one series transaction for both loads",
        ));
    }
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    let mut load = model.load.clone();
    let report = update(&mut load, cancelled)?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    Ok((load, report))
}

struct SeriesCandidate<R> {
    loads: [ApertureDgLoad; 2],
    reports: [R; 2],
    gas: GasState,
    impedance: f64,
}

fn prepare_series<R, C: FnMut() -> bool>(
    model: &ApertureDg<LoadedAperture>,
    mut update: impl FnMut(usize, &mut ApertureDgLoad, &mut C) -> Result<(R, f64), AcousticRealizeError>,
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut C,
) -> Result<SeriesCandidate<R>, AcousticRealizeError> {
    if model.supply_channel.is_some() {
        return Err(invalid("live DG series gas refresh must include the fixed supply passage's momentum and parameter work"));
    }
    let upstream = model
        .upstream
        .as_ref()
        .ok_or_else(|| invalid("live series gas refresh requires an original upstream DG tract"))?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    let mut loads = [upstream.clone(), model.load.clone()];
    let (up_report, up_work) = update(0, &mut loads[0], cancelled)?;
    let (down_report, down_work) = update(1, &mut loads[1], cancelled)?;
    let gas = *loads[1].duct().inlet_gas();
    let impedance = loads[0].impedance_pa_s_m3() + loads[1].impedance_pa_s_m3();
    let absolute_work = up_work + down_work;
    if loads[0].duct().inlet_gas() != &gas
        || !impedance.is_finite()
        || impedance <= 0.0
        || !maximum_absolute_parameter_work_j.is_finite()
        || maximum_absolute_parameter_work_j < 0.0
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
    Ok(SeriesCandidate {
        loads,
        reports: [up_report, down_report],
        gas,
        impedance,
    })
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
    ) -> Result<DgPortGasUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            |load, cancel| {
                load.rebind_gas_field_controlled(
                    field,
                    refresh,
                    maximum_absolute_parameter_work_j,
                    cancel,
                )
            },
            cancelled,
        )?;
        // This update validates before touching the original source. There is
        // no fallible operation between source and complete load publication.
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.load = load;
        Ok(report)
    }

    /// Refresh selected original bore/branch sources and the same moving valve.
    /// Branch-only changes update the exact mean load impedance while retaining
    /// the original inlet gas and source viscosity. Every mechanical/contact
    /// coordinate, material history, terminal state and sample clock survives.
    /// The native regional transaction owns all parameter-work accounting.
    pub fn rebind_regional_gas_fields_controlled(
        &mut self,
        main: Option<(&AxialGasField, DgGasRefresh)>,
        branches: &[DgBranchGasInput<'_>],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgRegionalGasUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            |load, cancel| {
                load.rebind_regional_gas_fields_controlled(
                    main,
                    branches,
                    maximum_absolute_parameter_work_j,
                    cancel,
                )
            },
            cancelled,
        )?;
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
    ) -> Result<DgPortGasUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            |load, cancel| {
                load.rebind_gas_field_controlled(
                    field,
                    refresh,
                    maximum_absolute_parameter_work_j,
                    cancel,
                )
            },
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(load.duct().inlet_gas(), load.impedance_pa_s_m3())?;
        self.load = load;
        Ok(report)
    }

    /// The same regional source transaction for original single, independent
    /// paired or material-modal mechanics. An attached tract requires the series
    /// regional API so both complete loads and their source publish together.
    pub fn rebind_regional_gas_fields_controlled(
        &mut self,
        main: Option<(&AxialGasField, DgGasRefresh)>,
        branches: &[DgBranchGasInput<'_>],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgRegionalGasUpdate, AcousticRealizeError> {
        let (load, report) = prepare_single(
            self,
            |load, cancel| {
                load.rebind_regional_gas_fields_controlled(
                    main,
                    branches,
                    maximum_absolute_parameter_work_j,
                    cancel,
                )
            },
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
    ) -> Result<[DgPortGasUpdate; 2], AcousticRealizeError> {
        let candidate = prepare_series(
            self,
            |i, load, cancel| {
                let report = load.rebind_gas_field_controlled(
                    fields[i],
                    refresh[i],
                    maximum_absolute_parameter_work_j,
                    cancel,
                )?;
                let work = report.absolute_parameter_work_j;
                Ok((report, work))
            },
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(&candidate.gas, candidate.impedance)?;
        let [up, down] = candidate.loads;
        self.upstream = Some(up);
        self.load = down;
        Ok(candidate.reports)
    }

    /// Refresh named original regions on either tract and their common source.
    /// Inputs and reports are upstream/downstream; an unselected side retains
    /// its exact owners and reports `None`. At least one source must be selected.
    /// Complete inlet gases must still match, and one allowance bounds the SUM
    /// of absolute native main/branch/channel work across both sides. Histories
    /// and mechanics are published only after the common-inlet and source checks.
    pub fn rebind_series_regional_gas_fields_controlled(
        &mut self,
        main: [Option<(&AxialGasField, DgGasRefresh)>; 2],
        branches: [&[DgBranchGasInput<'_>]; 2],
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Option<DgRegionalGasUpdate>; 2], AcousticRealizeError> {
        if main.iter().all(Option::is_none) && branches.iter().all(|inputs| inputs.is_empty()) {
            return Err(invalid(
                "regional tract refresh requires a selected original source",
            ));
        }
        let candidate = prepare_series(
            self,
            |i, load, cancel| {
                if main[i].is_none() && branches[i].is_empty() {
                    return Ok((None, 0.0));
                }
                let report = load.rebind_regional_gas_fields_controlled(
                    main[i],
                    branches[i],
                    maximum_absolute_parameter_work_j,
                    cancel,
                )?;
                let work = report.absolute_parameter_work_j;
                Ok((Some(report), work))
            },
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        self.aperture
            .rebind_exact_gas_load(&candidate.gas, candidate.impedance)?;
        let [up, down] = candidate.loads;
        self.upstream = Some(up);
        self.load = down;
        Ok(candidate.reports)
    }
}
