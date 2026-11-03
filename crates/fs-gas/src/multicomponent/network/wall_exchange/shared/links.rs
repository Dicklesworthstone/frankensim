//! Ordered thermal ports on evolving mixtures' original shared-wall states.

use super::*;
pub(super) use crate::multicomponent::wall_exchange::shared::links::WallLinks;
pub use crate::multicomponent::wall_exchange::shared::links::{
    SharedWallLink, SharedWallLinkReport,
};

impl NetworkWallOperator<'_> {
    /// Connect the bound original wall inventories through passive heat ports.
    ///
    /// Bind after `with_shared_walls`, before a source interval. Each source
    /// half first exchanges heat between the gas and its shared wall, then
    /// solves these wall-to-wall ports in authored order. Chemistry and flow
    /// between the halves change the gas response used by the next half.
    ///
    /// Each port solves `Q = dt G (Ta(-Q) - Tb(Q))` on the original h/u or
    /// capacity laws. Transfer is internal, without a second wall inventory or
    /// gas heat source. This is first-order composition, not a simultaneous
    /// thermal graph solve; timestep refinement must control link-order error.
    /// The existing source-query budget covers gas and link stages together.
    pub fn with_shared_wall_links(
        mut self,
        links: Vec<SharedWallLink>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        let shared = self.shared.as_mut().ok_or(WallExchangeError::Invalid(
            "wall links require shared inventories",
        ))?;
        if self.last_report.is_some()
            || self.trial_active
            || shared.links.is_some()
            || links.is_empty()
            || links.len() > shared.work.limit / 2
        {
            return Err(WallExchangeError::Invalid(
                "wall link binding state or work budget",
            ));
        }
        for link in &links {
            poll(&mut cancelled)?;
            if link.walls[0] == link.walls[1]
                || link.walls.iter().any(|&i| i >= shared.groups.len())
                || !link.conductance_w_per_k.is_finite()
                || link.conductance_w_per_k < 0.0
            {
                return Err(WallExchangeError::Invalid(
                    "wall link indices or conductance",
                ));
            }
        }
        // A standalone source has one stage; a complete network transaction
        // has two source halves. Retain both sets of actual link endpoints.
        let capacity = links
            .len()
            .checked_mul(2)
            .ok_or(WallExchangeError::Invalid("wall link report capacity"))?;
        let mut candidate = Vec::new();
        let mut accepted = Vec::new();
        candidate
            .try_reserve_exact(capacity)
            .map_err(|_| WallExchangeError::Invalid("wall link report allocation"))?;
        accepted
            .try_reserve_exact(capacity)
            .map_err(|_| WallExchangeError::Invalid("wall link report allocation"))?;
        shared.links = Some(WallLinks {
            ports: links,
            candidate,
            accepted,
            #[cfg(feature = "material")]
            conduction: None,
        });
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original authored port order and conductances, unchanged by evolution.
    pub fn shared_wall_links(&self) -> Option<&[SharedWallLink]> {
        Some(&self.shared.as_ref()?.links.as_ref()?.ports)
    }

    /// Every link stage in the last fully accepted source transaction.
    ///
    /// One chunk in authored link order is retained per source half: a direct
    /// `advance_for` has one chunk; a complete flow/source composition has two.
    /// Each receipt keeps its actual half duration and endpoint states. Later
    /// stages never overwrite earlier flux equations with the final wall state.
    /// Failure retains the entire preceding accepted report sequence.
    pub fn shared_wall_link_reports(&self) -> Option<&[SharedWallLinkReport]> {
        let shared = self.shared.as_ref()?;
        shared.last_report?;
        Some(&shared.links.as_ref()?.accepted)
    }

    /// Bind original bulk conductivity to either side of existing wall links.
    ///
    /// Rows follow link order, with each A/B material attached to that wall's
    /// original state and caloric source. On a bound link the authored G is
    /// the finite interface/contact conductance; selected L/(k A) resistances
    /// add once in series. Unbound links retain complete G, and zero G remains
    /// exact insulation. Perfect contact and spatial temperature profiles are
    /// not inferred by this finite-interface, uniform-inventory model.
    ///
    /// Initial-state layers refresh at the actual start of each link stage;
    /// implicit layers refresh per candidate without an invented derivative or
    /// an initial-flux cap. Later links must leave every earlier material source
    /// supported before publication. Both source halves share original storage,
    /// source-work accounting, cancellation and the network's atomic commit.
    #[cfg(feature = "material")]
    pub fn with_shared_wall_link_conduction(
        mut self,
        layers: Vec<[Option<MaterialWallConduction>; 2]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, WallExchangeError> {
        poll(&mut cancelled)?;
        if self.last_report.is_some() || self.trial_active {
            return Err(WallExchangeError::Invalid(
                "wall link material must bind before evolution",
            ));
        }
        let shared = self.shared.as_mut().ok_or(WallExchangeError::Invalid(
            "wall link material requires shared walls",
        ))?;
        let links = shared.links.as_mut().ok_or(WallExchangeError::Invalid(
            "wall link material requires existing links",
        ))?;
        links.bind_conduction(layers, |i| &shared.groups[i].storage, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(self)
    }

    /// Original material sources and resistance geometry in this link's A/B order.
    #[cfg(feature = "material")]
    pub fn shared_wall_link_conduction(
        &self,
        link: usize,
    ) -> Option<&[Option<MaterialWallConduction>; 2]> {
        self.shared
            .as_ref()?
            .links
            .as_ref()?
            .conduction_layers(link)
    }
}
