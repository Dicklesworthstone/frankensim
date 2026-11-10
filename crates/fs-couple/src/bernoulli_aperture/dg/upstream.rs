//! The owner of the source side upstream of a DG valve: an original DG tract, or a player's airway.
//!
//! [`UpstreamOwner::Dg`] is the original exact DG tract of [`super::tract`]: its inlet points away from the valve and
//! it owns its acoustic, branch, wall and terminal ledger. [`UpstreamOwner::Airway`] is the player's airway of
//! [`crate::render::player_tract`] (a measured area function realized as a passive wave-variable FIR,
//! [`crate::upstream_airway`]). Both give the valve the same law, `p = H - Z U` against the flow `U` into it, and
//! the valve is declared against `Z_down + Z_up` for either. They differ in who keeps the books:
//!
//! - the DG tract keeps its own stored energy, dissipation and terminal exchange, which the DG frame adds to its
//!   totals ([`super::DgTractFrame`]);
//! - the airway keeps NO stored energy in the DG account. Its FIRs are a lumped, band-limited port law, so it owns
//!   an absorbed-work ledger and a passivity certificate instead ([`PlayerTract::dg_commit`]), and the DG frame
//!   reports its delivered and absorbed work as a separate term ([`super::DgApertureFrame::airway`]), as it reports
//!   `supply`. [`UpstreamOwner::stored_energy_j`] is therefore zero for an airway, by definition and not by default.
//!
//! The constriction noise of a card is a series source at the valve for the airway (the `series_pa` of
//! [`UpstreamOwner::Airway`]), because the airway's far end takes the gesture pressure alone; for a DG tract it
//! joins the gesture pressure at the tract inlet, as before.

use super::{ApertureDgLoad, invalid, series_inlet_state};
use crate::acoustic_realize::AcousticRealizeError;
use crate::render::player_tract::PlayerTract;
use fs_material::gas::GasState;

/// The owner upstream of a DG valve (module doc).
#[derive(Clone)]
pub enum UpstreamOwner {
    /// An original exact DG tract whose inlet points away from the valve.
    Dg(ApertureDgLoad),
    /// A player's airway between the card's blowing pressure and the valve.
    Airway {
        /// The realized airway, its one-pole, its absorbed-work ledger and passivity certificate.
        tract: Box<PlayerTract>,
        /// The series source at the valve for the NEXT step, taken (reset to zero) by the step that accepts it [Pa].
        series_pa: f64,
    },
}

impl From<ApertureDgLoad> for UpstreamOwner {
    fn from(load: ApertureDgLoad) -> Self {
        Self::Dg(load)
    }
}

impl From<PlayerTract> for UpstreamOwner {
    fn from(tract: PlayerTract) -> Self {
        Self::Airway { tract: Box::new(tract), series_pa: 0.0 }
    }
}

/// The borrowed form a series step works on, so a private copy of a DG load (the finite-gas path) and the owner kept
/// in the valve take the same step code.
pub(super) enum UpstreamRef<'a> {
    Dg(&'a mut ApertureDgLoad),
    Airway { tract: &'a mut PlayerTract, series_pa: f64 },
}

impl UpstreamOwner {
    /// The original DG tract, if that is what sits upstream; `None` for an airway.
    #[must_use]
    pub const fn as_dg(&self) -> Option<&ApertureDgLoad> {
        match self {
            Self::Dg(load) => Some(load),
            Self::Airway { .. } => None,
        }
    }

    /// The player's airway, if that is what sits upstream; `None` for a DG tract.
    #[must_use]
    pub fn as_airway(&self) -> Option<&PlayerTract> {
        match self {
            Self::Dg(_) => None,
            Self::Airway { tract, .. } => Some(tract),
        }
    }

    pub(super) fn as_ref_mut(&mut self) -> UpstreamRef<'_> {
        match self {
            Self::Dg(load) => UpstreamRef::Dg(load),
            Self::Airway { tract, series_pa } => UpstreamRef::Airway { tract, series_pa: *series_pa },
        }
    }

    /// Stored energy in the DG account [J]: the DG tract's own; ZERO for an airway, whose FIR port law keeps its
    /// work in its ledger instead (module doc).
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        match self {
            Self::Dg(load) => load.stored_energy_j(),
            Self::Airway { .. } => 0.0,
        }
    }

    /// The series impedance this owner puts against the flow into the valve [Pa s/m^3].
    #[must_use]
    pub fn impedance_pa_s_m3(&self) -> f64 {
        match self {
            Self::Dg(load) => load.impedance_pa_s_m3(),
            Self::Airway { tract, .. } => tract.series_impedance_pa_s_m3(),
        }
    }

    /// Resolve the cold inlet gas and the SUM of the owner's and the bore's mean impedances: what the valve is
    /// declared against. A DG tract follows [`series_inlet_state`]. An airway must be cold, realized in the bore's
    /// complete inlet gas on the bore's own clock, and the bore must be cold with no driven exterior mouth; the
    /// impedance is the airway's series impedance plus the bore's, in that order, as for a DG tract.
    ///
    /// # Errors
    /// Any of those conditions failing, or a non-finite or non-positive summed impedance.
    pub fn inlet_state(
        &self,
        downstream: &ApertureDgLoad,
    ) -> Result<(GasState, f64), AcousticRealizeError> {
        match self {
            Self::Dg(load) => series_inlet_state(load, downstream),
            Self::Airway { tract, .. } => {
                let gas = *downstream.duct().inlet_gas();
                let impedance = tract.series_impedance_pa_s_m3() + downstream.impedance_pa_s_m3();
                if downstream.has_driven_mouth()
                    || downstream.state().step_count() != 0
                    || !tract.is_cold()
                    || tract.gas() != &gas
                    || (1.0 / f64::from(tract.airway().sample_rate_hz())).to_bits()
                        != downstream.propagator().dt().to_bits()
                    || !impedance.is_finite()
                    || impedance <= 0.0
                {
                    return Err(invalid(
                        "airway tract requires a cold airway realized in the bore's complete inlet gas on the bore's clock, a cold bore with no driven exterior mouth and a finite summed impedance",
                    ));
                }
                Ok((gas, impedance))
            }
        }
    }
}
