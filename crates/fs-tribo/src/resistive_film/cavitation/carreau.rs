//! Original thermal Carreau transport in the retained Elrod--Adams inventory.
//!
//! The original slit source resolves pressure and relative sliding together;
//! mean entrainment carries upwind liquid content. Both constitutive families
//! share the existing nonlinear pressure/content solve and independent mass
//! audit. Source temperature is frozen during a step and changes only through
//! an explicit source replacement.
//!
//! This model returns retained liquid mass, face flow and endpoint NORMAL
//! pressure force. A full-gap source profile does not determine the geometry
//! or heat of liquid streamers in a partially filled film. Cavitated sliding
//! traction, viscous heat partition and thermal evolution therefore remain
//! outside this owner; no full-gap heat is relabeled as cavitated heat.

use fs_material::fluid_transport::rheology::slit::{CarreauSlitControls, LaminarCarreauSlit};

use super::{
    CavitatingFilm, CavitatingFilmState, CavitationLiquid, CavitationReservoir,
    CavitationSolveConfig, CavitationStep, FilmError, ResistiveFilm, dot, fail, poll,
};
use crate::resistive_film::carreau::{CarreauFilm, CarreauFilmControls, WallVelocityPorts};

/// Original slit-work controls and bounded pressure/content backtracking.
/// Iteration count and mass tolerances belong to `CavitationSolveConfig`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CarreauCavitationControls {
    pub passage: CarreauSlitControls,
    pub maximum_line_search_steps: usize,
}

impl Default for CarreauCavitationControls {
    fn default() -> Self {
        Self {
            passage: CarreauSlitControls::default(),
            maximum_line_search_steps: 48,
        }
    }
}

/// Mass-conserving cavitation with the original thermal Carreau material law.
///
/// Each passage retains its own prepared temperature, transition, plateaus and
/// supported shear-rate ceiling. Width and length must match the original
/// graph; the slit aperture is rebound to the actual affine channel gap.
/// Pressure-dependent viscosity and compressibility of trapped full liquid
/// are not added. All states retain the native inventory identity and mass
/// audit; the liquid content is not a vapor-mass or latent-heat coordinate.
#[derive(Clone, Debug)]
pub struct CavitatingCarreauFilm {
    inventory: CavitatingFilm,
    rheology: CarreauFilm,
}

/// Matching endpoint source owner and conservative candidate, accepted together.
#[derive(Clone, Debug)]
pub struct CavitatingCarreauContinuation {
    pub film: CavitatingCarreauFilm,
    pub step: CavitationStep,
}

impl CavitatingFilm {
    /// Bind original Carreau passages without resetting accepted liquid mass.
    /// Geometry, explicit supplies, mass controls and inventory identity survive.
    /// A native Barus law refuses through the original source consumer.
    pub fn with_carreau(
        &self,
        passages: Vec<LaminarCarreauSlit>,
        controls: CarreauCavitationControls,
    ) -> Result<CavitatingCarreauFilm, FilmError> {
        let rheology = prepare_rheology(self.geometry(), passages, controls)?;
        let inventory =
            self.with_channel_viscosities(rheology.pressure_operator().channel_viscosities_pa_s())?;
        Ok(CavitatingCarreauFilm {
            inventory,
            rheology,
        })
    }
}

impl CavitatingCarreauFilm {
    /// Bind source passages to the existing native liquid-content owner.
    pub fn new(
        geometry: ResistiveFilm,
        passages: Vec<LaminarCarreauSlit>,
        liquid: CavitationLiquid,
        solve: CavitationSolveConfig,
        controls: CarreauCavitationControls,
    ) -> Result<Self, FilmError> {
        let rheology = prepare_rheology(&geometry, passages, controls)?;
        let inventory = CavitatingFilm::new(rheology.pressure_operator().clone(), liquid, solve)?;
        Ok(Self {
            inventory,
            rheology,
        })
    }

    /// Original graph with prepared zero-shear coefficients. Evaluating this
    /// geometry alone does not represent nonlinear or cavitating transport.
    pub fn geometry(&self) -> &ResistiveFilm {
        self.inventory.geometry()
    }

    pub fn liquid(&self) -> CavitationLiquid {
        self.inventory.liquid()
    }

    pub fn passages(&self) -> &[LaminarCarreauSlit] {
        self.rheology.passages()
    }

    pub fn controls(&self) -> CarreauCavitationControls {
        CarreauCavitationControls {
            passage: self.rheology.controls().passage,
            maximum_line_search_steps: self.rheology.controls().pressure.maximum_line_search_steps,
        }
    }

    pub fn reservoir_conditions(
        &self,
    ) -> impl ExactSizeIterator<Item = Option<CavitationReservoir>> + '_ {
        self.inventory.reservoir_conditions()
    }

    /// Explicit absolute-pressure and incoming-content supplies, preserving mass.
    pub fn with_reservoirs(
        &self,
        reservoirs: &[Option<CavitationReservoir>],
    ) -> Result<Self, FilmError> {
        Ok(Self {
            inventory: self.inventory.with_reservoirs(reservoirs)?,
            rheology: self.rheology.clone(),
        })
    }

    /// Install explicit endpoint temperatures/laws with compatible inventory.
    /// This does not advance time, temperature, density or rupture pressure.
    pub fn with_passages(&self, passages: Vec<LaminarCarreauSlit>) -> Result<Self, FilmError> {
        let rheology = self.rheology.with_passages(passages)?;
        let inventory = self
            .inventory
            .with_channel_viscosities(rheology.pressure_operator().channel_viscosities_pa_s())?;
        Ok(Self {
            inventory,
            rheology,
        })
    }

    pub fn initial_state(&self, q: &[f64]) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initial_state(q)
    }

    /// Explicit starved content, at rupture pressure when partially filled.
    pub fn initialize(
        &self,
        q: &[f64],
        liquid_fraction: &[f64],
    ) -> Result<CavitatingFilmState, FilmError> {
        self.inventory.initialize(q, liquid_fraction)
    }

    /// Advance mass with the actual simultaneous pressure/sliding profile.
    ///
    /// `entrainment_m3_s` is the signed full-film mean-wall volume flow, which
    /// the native owner multiplies by upwind content. The separately authored
    /// `relative_wall_speed_m_s` is upper minus lower speed in that SAME channel
    /// orientation. No-slip walls give entrainment `width*gap*(Ulower+Uupper)/2`.
    /// The caller owns wall-time quadrature and thin-gap creeping-flow validity.
    ///
    /// Constant sources also use the shared nonlinear engine so original rate
    /// limits are checked at every accepted face. All failed source, budget,
    /// pressure, mass or cancellation checks leave accepted state untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        poll(&mut cancelled)?;
        self.validate_relative_motion(q_next, relative_wall_speed_m_s)?;
        self.advance_kernel(
            &self.inventory,
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            &mut cancelled,
        )
    }

    /// Continue OLD liquid mass using explicit endpoint liquid and source data.
    ///
    /// The native balance remains
    /// `rho_next*(V_next*theta_next + dt*div Q_next) = rho_old*V_old*theta_old`.
    /// A changed density never reprices the old inventory. Explicit supplies
    /// must remain valid at the new rupture pressure. A liquid-data change
    /// creates a distinct inventory identity; a source-only change preserves it.
    ///
    /// This first-order endpoint constitutive split supplies no caloric or
    /// latent energy. The surrounding owner accepts the returned source and
    /// candidate together and owns temperature, material history and refinement.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_liquid(
        &self,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        liquid: CavitationLiquid,
        passages: Vec<LaminarCarreauSlit>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CavitatingCarreauContinuation, FilmError> {
        poll(&mut cancelled)?;
        self.inventory.validate_state(state)?;
        self.validate_relative_motion(q_next, relative_wall_speed_m_s)?;
        let rheology = self.rheology.with_passages(passages)?;
        let inventory = self.inventory.prepare_liquid(
            liquid,
            rheology.pressure_operator().channel_viscosities_pa_s(),
        )?;
        let film = Self {
            inventory,
            rheology,
        };
        let step = film.advance_kernel(
            &self.inventory,
            state,
            q_next,
            entrainment_m3_s,
            relative_wall_speed_m_s,
            dt_s,
            &mut cancelled,
        )?;
        Ok(CavitatingCarreauContinuation { film, step })
    }

    fn validate_relative_motion(&self, q: &[f64], relative: &[f64]) -> Result<(), FilmError> {
        if relative.len() != self.inventory.channel_count()
            || relative.iter().any(|v| !v.is_finite())
        {
            return fail("invalid cavitating Carreau relative wall speeds");
        }
        self.geometry().validate_configuration(q)?;
        for (edge, &speed) in self.geometry().channels.iter().zip(relative) {
            if edge.gap.reference_m - dot(&edge.gap.closure, q) <= 0.0 && speed != 0.0 {
                return fail("relative sliding through a closed cavitating passage");
            }
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_kernel(
        &self,
        previous: &CavitatingFilm,
        state: &CavitatingFilmState,
        q_next: &[f64],
        entrainment_m3_s: &[f64],
        relative_wall_speed_m_s: &[f64],
        dt_s: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<CavitationStep, FilmError> {
        self.inventory.advance_kernel(
            previous,
            state,
            q_next,
            entrainment_m3_s,
            dt_s,
            Some(super::power_law::Transport::Carreau {
                rheology: &self.rheology,
                relative: relative_wall_speed_m_s,
            }),
            cancelled,
        )
    }
}

fn prepare_rheology(
    geometry: &ResistiveFilm,
    passages: Vec<LaminarCarreauSlit>,
    controls: CarreauCavitationControls,
) -> Result<CarreauFilm, FilmError> {
    // This owner authors wall speeds and reservoir contents at each mass step.
    // Bind the shared original material consumer without importing unrelated
    // full-film wall rows or gauge-pressure supplies into those native ports.
    let reservoirs: Vec<_> = geometry
        .channels
        .iter()
        .map(|edge| edge.to.is_none().then_some(0.0))
        .collect();
    let walls = vec![
        WallVelocityPorts {
            lower: vec![0.0; geometry.port_count()],
            upper: vec![0.0; geometry.port_count()],
        };
        geometry.channel_count()
    ];
    let mut source_controls = CarreauFilmControls::default();
    source_controls.passage = controls.passage;
    source_controls.pressure.maximum_line_search_steps = controls.maximum_line_search_steps;
    CarreauFilm::new(geometry, passages, &reservoirs, walls, source_controls)
}
