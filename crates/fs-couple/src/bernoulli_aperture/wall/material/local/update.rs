//! Source-backed wall changes on a continuing native DG trajectory.

use super::*;
use fs_duct::dg::{DgDuct, DgPropagator, DgState, DgWallMaterialUpdate};

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

fn prepare(
    duct: &DgDuct,
    state: &DgState,
    propagator: &DgPropagator,
    walls: &[LocalMaterialWall],
    retained_band_hz: [f64; 2],
    maximum_absolute_parameter_work_j: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedUpdate, AcousticRealizeError> {
    if cancelled() {
        return Err(AcousticRealizeError::Cancelled);
    }
    if walls.is_empty() || walls.len() != duct.mesh().elements {
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
    for (element, wall) in walls.iter().enumerate() {
        if cancelled() {
            return Err(AcousticRealizeError::Cancelled);
        }
        if matches!(&wall.loss, WallMaterialLoss::Maxwell(_)) {
            return Err(invalid(
                "live hereditary wall updates need a source-identified branch-history transfer",
            ));
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
    let account = next_duct
        .rebind_compliant_wall_regions(
            &mut next_state,
            &pins,
            maximum_absolute_parameter_work_j,
            cancelled,
        )
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

impl MaterialWallDg {
    /// Replace a live uniform elastic/Kelvin--Voigt wall from original sources.
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
    /// frequency band. Native geometry, work-budget and Maxwell-history refusals
    /// apply. This intervention does not solve heat, expansion, changing gas,
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
        let proposal = prepare(
            &self.duct,
            state,
            propagator,
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
    /// Replace each original element's live material without averaging sources.
    ///
    /// This is the regional counterpart of `MaterialWallDg::rebind_controlled`:
    /// fixed geometry/gas, preserved wall displacement/momentum and gas-memory
    /// history, one separately budgeted parameter-work update, and atomic source,
    /// state and propagator publication. Replacement source bands must all cover
    /// the existing common band. Live changes to retained Maxwell histories
    /// refuse; original cold bindings may retain independent regional spectra.
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
