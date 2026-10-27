//! Source-backed wall changes on a continuing native DG trajectory.

use super::*;
use fs_duct::dg::{
    DgDuct, DgGasMaterialUpdate, DgPropagator, DgState, DgWallMaterialUpdate,
    DgWallRelaxationTransfer,
};
use fs_material::gas_field::AxialGasField;

fn duct_error(error: fs_duct::DuctError) -> AcousticRealizeError {
    match error {
        fs_duct::DuctError::Cancelled { .. } => AcousticRealizeError::Cancelled,
        error => AcousticRealizeError::Nonlinear(error.to_string()),
    }
}

struct PreparedUpdate {
    duct: DgDuct,
    state: DgState,
    propagator: DgPropagator,
    account: DgWallMaterialUpdate,
    band_hz: [f64; 2],
}

#[allow(clippy::too_many_arguments)]
fn prepare(
    duct: &DgDuct,
    state: &DgState,
    propagator: &DgPropagator,
    previous_walls: &[LocalMaterialWall],
    walls: &[LocalMaterialWall],
    retained_band_hz: [f64; 2],
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedUpdate, AcousticRealizeError> {
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    if walls.is_empty()
        || walls.len() != duct.mesh().elements
        || previous_walls.len() != walls.len()
    {
        return Err(invalid(
            "live DG wall update needs one source per original element",
        ));
    }
    let radius = duct.uniform_circular_radius_m().ok_or_else(|| {
        invalid("live material wall update requires the original circular geometry")
    })?;
    let nodes = duct.node_positions();
    let local = duct.mesh().order + 1;
    let mut band_hz = walls[0].band_hz;
    let mut pins = Vec::with_capacity(walls.len());
    let mut correspondence = Vec::new();
    let mut has_relaxation = false;
    for (element, wall) in walls.iter().enumerate() {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        let previous = &previous_walls[element];
        match (&previous.loss, &wall.loss) {
            (WallMaterialLoss::Maxwell(old), WallMaterialLoss::Maxwell(new)) => {
                has_relaxation = true;
                if previous.material.card_identity() != wall.material.card_identity()
                    || previous.thickness_m.to_bits() != wall.thickness_m.to_bits()
                    || previous.axial != wall.axial
                {
                    return Err(invalid(
                        "live Maxwell history requires the original material card, thickness and restraint",
                    ));
                }
                for name in [
                    DENSITY_PROPERTY,
                    YOUNG_MODULUS_PROPERTY,
                    POISSON_RATIO_PROPERTY,
                ] {
                    let selected = |source: &LocalMaterialWall| {
                        source
                            .material
                            .property(name)
                            .map(|property| property.answer().receipt.selected)
                    };
                    if selected(previous) != selected(wall) {
                        return Err(invalid(
                            "live Maxwell history must retain the original base material claims",
                        ));
                    }
                }
                correspondence.extend(old.history_correspondence(new)?.into_iter().map(
                    |(previous_arm, replacement_arm)| DgWallRelaxationTransfer {
                        element,
                        previous_arm,
                        replacement_arm,
                    },
                ));
            }
            (WallMaterialLoss::Maxwell(_), _) | (_, WallMaterialLoss::Maxwell(_)) => {
                return Err(invalid(
                    "live Maxwell update cannot introduce or discard hereditary source histories",
                ));
            }
            _ => {}
        }
        band_hz[0] = band_hz[0].max(wall.band_hz[0]);
        band_hz[1] = band_hz[1].min(wall.band_hz[1]);
        let base = element * local;
        pins.push(
            wall.base_patch(radius, nodes[base + local - 1] - nodes[base])?
                .wall,
        );
    }
    if band_hz[0] > retained_band_hz[0] || band_hz[1] < retained_band_hz[1] {
        return Err(invalid(
            "replacement wall sources must cover the existing live applicability band",
        ));
    }
    let mut next_state = state.clone();
    // The supplied stepping owner must actually belong to this accepted state.
    propagator.sync_nodal(&mut next_state).map_err(duct_error)?;
    let mut next_duct = duct.clone();
    let account = if has_relaxation {
        let spectra: Vec<_> = walls
            .iter()
            .map(|wall| match &wall.loss {
                WallMaterialLoss::Maxwell(source) => Some(source.material()),
                _ => None,
            })
            .collect();
        next_duct.rebind_compliant_wall_relaxation_regions(
            &mut next_state,
            &pins,
            &spectra,
            &correspondence,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    } else {
        next_duct.rebind_compliant_wall_regions(
            &mut next_state,
            &pins,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }
    .map_err(duct_error)?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    // The original certified exponential owner admits the new coefficients at
    // the SAME physical clock. A failed certificate publishes no source/state.
    let next_propagator = next_duct
        .propagator(propagator.dt(), propagator.tolerance())
        .map_err(duct_error)?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    Ok(PreparedUpdate {
        duct: next_duct,
        state: next_state,
        propagator: next_propagator,
        account,
        band_hz,
    })
}

/// Prepare the same live gas transfer and certified clock for either wall
/// source layout. Native wall coordinates and all immutable source assignments
/// remain with their existing owner; the new field retains its gas catalog.
fn prepare_gas(
    duct: &DgDuct,
    state: &DgState,
    propagator: &DgPropagator,
    field: &AxialGasField,
    viscothermal_error: Option<f64>,
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(DgDuct, DgState, DgPropagator, DgGasMaterialUpdate), AcousticRealizeError> {
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    let mut next_state = state.clone();
    propagator.sync_nodal(&mut next_state).map_err(duct_error)?;
    let mut next_duct = duct.clone();
    let report = if let Some(error) = viscothermal_error {
        next_duct.rebind_uniform_viscothermal_gas_field(
            &mut next_state,
            field,
            error,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    } else {
        next_duct.rebind_uniform_gas_field(
            &mut next_state,
            field,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
    }.map_err(duct_error)?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    let next_propagator = next_duct
        .propagator(propagator.dt(), propagator.tolerance())
        .map_err(duct_error)?;
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    Ok((next_duct, next_state, next_propagator, report))
}

impl MaterialWallDg {
    /// Refresh the continuing uniform gas from its complete shared source field.
    ///
    /// Native gas compression and momentum are preserved; pressure, volume flow,
    /// boundary/side-port efforts and reciprocal wall loading follow the new
    /// source T/pressure/composition. Wall source, displacement, velocity and
    /// every Maxwell history survive. The exact physical storage change enters
    /// the separate parameter-work ledger once, without advancing the clock.
    /// The owner changing the background must settle this signed work; the
    /// intervention does not solve heat, mass transport or mean-pressure traction.
    ///
    /// Field, native operator, state and certified propagator publish together.
    /// The original timestep/tolerance and wall applicability band are retained.
    /// Nonuniform gas and viscothermal GAS memory refuse, never average/reset;
    /// regional Maxwell wall memory remains supported by the native transfer.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (duct, next, stepping, report) = prepare_gas(
            &self.duct,
            state,
            propagator,
            field,
            None,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_duct = core::mem::replace(&mut self.duct, duct);
        let _old_state = core::mem::replace(state, next);
        let _old_propagator = core::mem::replace(propagator, stepping);
        Ok(report)
    }

    /// Refresh uniform source gas while retaining native circular diffusion.
    ///
    /// The original Bessel sections preserve every viscous arm's momentum and
    /// every thermal arm's compression volume. Wave and memory storage changes
    /// enter parameter work once; Maxwell wall sources and histories survive.
    /// The new gas must resolve the retained band's complex response and real
    /// loss under maximum_relative_error, without changing radius/pole support.
    /// The complete source field, native state and SAME-clock certified
    /// propagator publish together after all admission and cancellation checks.
    /// This is a prescribed material intervention, not a heat/species solver;
    /// its external owner owes the signed work and any background gas balance.
    pub fn rebind_viscothermal_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (duct, next, stepping, report) = prepare_gas(
            &self.duct,
            state,
            propagator,
            field,
            Some(maximum_relative_error),
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_duct = core::mem::replace(&mut self.duct, duct);
        let _old_state = core::mem::replace(state, next);
        let _old_propagator = core::mem::replace(propagator, stepping);
        Ok(report)
    }

    /// Replace a live uniform wall from original material sources.
    ///
    /// Supply a newly resolved material state (for example at a new supported
    /// temperature), thickness and restraint through `LocalMaterialWall`. The
    /// existing native transfer preserves radial displacement and conjugate
    /// momentum, and retains the complete acoustic/viscothermal history. Its
    /// signed storage change is charged once as parameter work; the external
    /// material/thermal owner owes that work and must not also call it heat.
    ///
    /// Source receipts, duct, accepted state and certified propagator publish
    /// together after every check. The propagator keeps its timestep/tolerance;
    /// no sample is consumed. A replacement must cover the existing declared
    /// frequency band. A Maxwell wall also preserves every arm's elastic
    /// extension and viscous offset, identified by its original property pair
    /// and selected claims on the SAME immutable card. Reordering source arms
    /// is supported; changing thickness/restraint, unlabelled authored spectra,
    /// new/dropped arms and zero/positive transitions refuse. A relaxation-time
    /// change alters future loss without instantaneous parameter work.
    /// This intervention does not solve heat, expansion, changing gas,
    /// or mass transfer. Cancellation is checked around the existing bounded
    /// dense propagator preparation, which supplies no internal cancel callback.
    pub fn rebind_controlled(
        &mut self,
        wall: LocalMaterialWall,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgWallMaterialUpdate, AcousticRealizeError> {
        let walls = vec![wall.clone(); self.duct.mesh().elements];
        let previous_walls = vec![self.wall.clone(); self.duct.mesh().elements];
        let proposal = prepare(
            &self.duct,
            state,
            propagator,
            &previous_walls,
            &walls,
            self.wall.band_hz,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_bound = core::mem::replace(
            self,
            Self {
                wall,
                duct: proposal.duct,
            },
        );
        let _old_state = core::mem::replace(state, proposal.state);
        let _old_propagator = core::mem::replace(propagator, proposal.propagator);
        Ok(proposal.account)
    }
}

impl RegionalMaterialWallDg {
    /// Source-retaining uniform gas refresh for independently sourced walls.
    /// Uses `MaterialWallDg::rebind_gas_field_controlled`'s physical transfer and
    /// atomic publication, retaining every original element's material receipt,
    /// wall motion and hereditary history without a regional average.
    pub fn rebind_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (duct, next, stepping, report) = prepare_gas(
            &self.duct,
            state,
            propagator,
            field,
            None,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_duct = core::mem::replace(&mut self.duct, duct);
        let _old_state = core::mem::replace(state, next);
        let _old_propagator = core::mem::replace(propagator, stepping);
        Ok(report)
    }

    /// Regional-wall counterpart of the uniform source wrapper's live
    /// viscothermal gas refresh. Every original wall receipt, local Maxwell
    /// history and applicability band is retained without regional averaging.
    pub fn rebind_viscothermal_gas_field_controlled(
        &mut self,
        field: &AxialGasField,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgGasMaterialUpdate, AcousticRealizeError> {
        let (duct, next, stepping, report) = prepare_gas(
            &self.duct,
            state,
            propagator,
            field,
            Some(maximum_relative_error),
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_duct = core::mem::replace(&mut self.duct, duct);
        let _old_state = core::mem::replace(state, next);
        let _old_propagator = core::mem::replace(propagator, stepping);
        Ok(report)
    }

    /// Replace each original element's live material without averaging sources.
    ///
    /// This is the regional counterpart of `MaterialWallDg::rebind_controlled`:
    /// fixed geometry/gas, preserved wall displacement/momentum and gas-memory
    /// history, one separately budgeted parameter-work update, and atomic source,
    /// state and propagator publication. Replacement source bands must all cover
    /// the existing common band. Maxwell histories keep the original local
    /// source-arm identity and support under the uniform wrapper's rules;
    /// homogeneous and regional storage layouts transfer without projection.
    pub fn rebind_controlled(
        &mut self,
        walls: Vec<LocalMaterialWall>,
        state: &mut DgState,
        propagator: &mut DgPropagator,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgWallMaterialUpdate, AcousticRealizeError> {
        let proposal = prepare(
            &self.duct,
            state,
            propagator,
            &self.walls,
            &walls,
            self.band_hz,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?;
        let _old_bound = core::mem::replace(
            self,
            Self {
                walls,
                band_hz: proposal.band_hz,
                duct: proposal.duct,
            },
        );
        let _old_state = core::mem::replace(state, proposal.state);
        let _old_propagator = core::mem::replace(propagator, proposal.propagator);
        Ok(proposal.account)
    }
}
