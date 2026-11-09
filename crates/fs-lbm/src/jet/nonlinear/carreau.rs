//! Original thermal Carreau-Yasuda viscosity in the native slender-jet step.
//!
//! In incompressible uniaxial extension, D=diag(-e/2,-e/2,e), e=du/dz.
//! Thus gamma=sqrt(2 D:D)=sqrt(3)|e|, tau_zz-tau_rr=3 mu(gamma)e,
//! and tau:D=3 mu(gamma)e². Axial force, dissipation and the viscous step
//! bound therefore consume one original material sample at each RK stage.
//!
//! This is the generalized-Newtonian extension of the original slender
//! equations, not a separately calibrated extensional viscosity or an elastic
//! normal stress. See the generalized-Newtonian discussion in Eggers &
//! Villermaux (2008), Physics of liquid jets, equation (320):
//! <https://people.maths.bris.ac.uk/~majge/rpp8_3_036601.pdf>.
//! The existing radius/slope/resolution limits remain essential near a neck.

use super::*;
use fs_material::fluid_transport::rheology::pipe::PreparedCarreau;
use fs_material::fluid_transport::rheology::{ThermalCarreauYasuda, ThermalRheologyState};

/// Fixed physical conditions and the original generalized-Newtonian law.
/// The temperature shift and calibrated rate support are never refitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauSlenderJetParameters {
    pub density_kg_m3: f64,
    pub surface_tension_n_m: f64,
    pub temperature_k: f64,
    pub rheology: ThermalCarreauYasuda,
}

/// The original calibration, optionally prepared at prescribed CELL temperatures.
/// Cell i owns the extensional rate between faces i-1/2 and i+1/2.
#[derive(Clone, Debug)]
pub(super) struct CarreauProfile {
    parameters: CarreauSlenderJetParameters,
    cells: Option<Vec<PreparedCarreau>>,
    uniform_newtonian: bool,
}

impl CarreauProfile {
    pub(super) fn uniform(parameters: CarreauSlenderJetParameters) -> Self {
        Self {
            uniform_newtonian: parameters.is_newtonian(),
            parameters,
            cells: None,
        }
    }

    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn try_clone_prepared(&self, gate: &CancelGate) -> Result<Self, NonlinearJetError> {
        let original = self.cells.as_ref().ok_or(NonlinearJetError::InvalidState(
            "thermal jet requires prepared cell materials",
        ))?;
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(original.len())
            .map_err(|_| NonlinearJetError::Allocation)?;
        for (i, &cell) in original.iter().enumerate() {
            poll_at(gate, i)?;
            cells.push(cell);
        }
        Ok(Self {
            parameters: self.parameters,
            cells: Some(cells),
            uniform_newtonian: self.uniform_newtonian,
        })
    }

    #[cfg(feature = "enthalpy-coupling")]
    pub(super) fn reprepare_temperature(
        &mut self,
        cell: usize,
        temperature_k: f64,
    ) -> Result<(), NonlinearJetError> {
        let cells = self.cells.as_mut().ok_or(NonlinearJetError::InvalidState(
            "thermal jet requires prepared cell materials",
        ))?;
        cells[cell] = PreparedCarreau::new(
            self.parameters.rheology,
            temperature_k,
            cells[cell].maximum_shear_rate_per_s(),
        )
        .map_err(|source| NonlinearJetError::Rheology {
            cell: Some(cell),
            source,
        })?;
        // Evolving temperatures need not remain uniform even for a linear
        // rate law. The coupled owner uses only explicit original stresses.
        self.uniform_newtonian = false;
        Ok(())
    }

    pub(super) const fn is_uniform_newtonian(&self) -> bool {
        self.uniform_newtonian
    }

    pub(super) const fn density_kg_m3(&self) -> f64 {
        self.parameters.density_kg_m3
    }

    pub(super) fn temperature_k(&self, cell: usize) -> f64 {
        self.cells
            .as_ref()
            .map_or(self.parameters.temperature_k, |cells| {
                cells[cell].temperature_k()
            })
    }

    pub(super) fn coefficients(
        &self,
        axial_rate_per_s: f64,
        cell: usize,
    ) -> Result<[f64; 2], NonlinearJetError> {
        if let Some(cells) = &self.cells {
            let rate = finite(3.0_f64.sqrt() * axial_rate_per_s.abs())?;
            cells[cell]
                .viscosity_and_tangent(rate)
                .map_err(|source| NonlinearJetError::Rheology {
                    cell: Some(cell),
                    source,
                })
        } else {
            let sample = self.parameters.sample(axial_rate_per_s, Some(cell))?;
            Ok([
                sample.dynamic_viscosity_pa_s(),
                sample.stress_tangent_pa_s(),
            ])
        }
    }

    fn sample(
        &self,
        axial_rate_per_s: f64,
        cell: usize,
    ) -> Result<ThermalRheologyState, NonlinearJetError> {
        if self.cells.is_none() {
            return self.parameters.sample(axial_rate_per_s, Some(cell));
        }
        // Admit the same physical envelope as the actual solver, then return
        // the complete original source response (including temperature slopes).
        self.coefficients(axial_rate_per_s, cell)?;
        let rate = finite(3.0_f64.sqrt() * axial_rate_per_s.abs())?;
        self.parameters
            .rheology
            .evaluate(self.temperature_k(cell), rate)
            .map_err(|source| NonlinearJetError::Rheology {
                cell: Some(cell),
                source,
            })
    }
}

impl CarreauSlenderJetParameters {
    pub(super) fn is_newtonian(&self) -> bool {
        let [mu0, mu_inf, time, index, _] = self.rheology.curve().parameters();
        time == 0.0 || index == 1.0 || mu0 == mu_inf
    }

    pub(super) fn sample(
        &self,
        axial_rate_per_s: f64,
        cell: Option<usize>,
    ) -> Result<ThermalRheologyState, NonlinearJetError> {
        let rate = finite(3.0_f64.sqrt() * axial_rate_per_s.abs())?;
        self.rheology
            .evaluate(self.temperature_k, rate)
            .map_err(|source| NonlinearJetError::Rheology { cell, source })
    }
}

impl PeriodicSlenderJet {
    /// Construct the same periodic jet with source-resolved Carreau viscosity.
    ///
    /// Every stress and loss evaluates the original law at the local
    /// sqrt(3)*|du/dz|, in inverse seconds. The source must support zero rate
    /// and every actual RK stage at this fixed temperature. Unsupported
    /// algebraic Newton line-search candidates may backtrack within the work
    /// budget; source refusal of an actual stage rolls the entire advance back.
    /// Neither path clips a rate or substitutes another viscosity. Density and surface tension are
    /// authored constants; this owner adds no caloric evolution, gravity,
    /// viscoelastic stress, gas drag or continuation through pinch-off.
    ///
    /// The native conservative flux supports explicit SSPRK3 and the original
    /// third-order IMEX viscosity policy. Nonlinear implicit stages retain the
    /// original material stress and differential tangent, with bounded Newton
    /// work selected by [`Self::set_implicit_viscosity_controls`]. Scratch is
    /// prepared on policy selection and reused by successful steps.
    pub fn new_carreau(
        parameters: CarreauSlenderJetParameters,
        length_m: f64,
        radii_m: &[f64],
        velocities_m_s: &[f64],
        limits: SlenderJetLimits,
        gate: &CancelGate,
    ) -> Result<Self, NonlinearJetError> {
        poll(gate)?;
        if !parameters.density_kg_m3.is_finite()
            || parameters.density_kg_m3 <= 0.0
            || !parameters.surface_tension_n_m.is_finite()
            || parameters.surface_tension_n_m < 0.0
        {
            return Err(NonlinearJetError::InvalidInput(
                "Carreau liquid coefficients",
            ));
        }
        let reference = parameters.sample(0.0, None)?;
        let native = SlenderJetParameters {
            density_kg_m3: parameters.density_kg_m3,
            surface_tension_n_m: parameters.surface_tension_n_m,
            kinematic_viscosity_m2_s: positive(
                reference.dynamic_viscosity_pa_s() / parameters.density_kg_m3,
            )?,
        };
        Self::new_constitutive(
            native,
            Some(CarreauProfile::uniform(parameters)),
            length_m,
            radii_m,
            velocities_m_s,
            limits,
            gate,
        )
    }

    /// Original calibration and initial uniform SI conditions, or None for the
    /// Newtonian path. A prescribed cell-temperature field overrides only its
    /// temperature; see [`Self::prepared_carreau_cells`] for actual held values.
    #[must_use]
    pub const fn carreau_parameters(&self) -> Option<&CarreauSlenderJetParameters> {
        match &self.carreau {
            Some(profile) => Some(&profile.parameters),
            None => None,
        }
    }

    /// Original source preparations at the prescribed cell temperatures.
    /// None means the original uniform temperature is still in use.
    #[must_use]
    pub fn prepared_carreau_cells(&self) -> Option<&[PreparedCarreau]> {
        self.carreau.as_ref()?.cells.as_deref()
    }

    /// Prescribe one fluid temperature at each original axial CELL.
    ///
    /// Cell i contains the extensional stress between its two adjacent velocity
    /// faces. Temperature is held throughout subsequent advance calls, until
    /// explicitly replaced here; it is not advected or evolved by this owner.
    /// Density and surface tension remain the original constants. This is a
    /// thermoviscous momentum solve, not thermal expansion or Marangoni flow.
    ///
    /// The original shift is prepared once per cell. The caller declares a
    /// common finite physical invariant-rate ceiling; every temperature must
    /// admit its full interval from rest through that ceiling. Neither rates
    /// nor original reference support are clipped. The ceiling is a constitutive
    /// work/domain limit, not a laminar-transition or accuracy certificate.
    ///
    /// The actual accepted strain field is checked before publication. Source,
    /// allocation or cancellation refusal preserves the old source field,
    /// numerical policy, physical state, diagnostics and time. Reusable solve
    /// scratch may be prepared. No accepted flow or source is reconstructed.
    pub fn set_carreau_cell_temperatures(
        &mut self,
        temperatures_k: &[f64],
        maximum_invariant_rate_per_s: f64,
        gate: &CancelGate,
    ) -> Result<(), NonlinearJetError> {
        poll(gate)?;
        let parameters = self
            .carreau
            .as_ref()
            .ok_or(NonlinearJetError::InvalidInput(
                "cell temperatures require an original Carreau jet",
            ))?
            .parameters;
        let n = self.accepted.a.len();
        if temperatures_k.len() != n {
            return Err(NonlinearJetError::InvalidInput(
                "jet cell temperature count",
            ));
        }
        let mut cells = Vec::new();
        cells
            .try_reserve_exact(n)
            .map_err(|_| NonlinearJetError::Allocation)?;
        for (i, &temperature) in temperatures_k.iter().enumerate() {
            poll_at(gate, i)?;
            let prepared = PreparedCarreau::new(
                parameters.rheology,
                temperature,
                maximum_invariant_rate_per_s,
            )
            .map_err(|source| NonlinearJetError::Rheology {
                cell: Some(i),
                source,
            })?;
            cells.push(prepared);
        }
        let reference = positive(cells[0].zero_shear_viscosity_pa_s() / parameters.density_kg_m3)?;
        let uniform_newtonian = parameters.is_newtonian()
            && cells.iter().all(|cell| {
                cell.zero_shear_viscosity_pa_s() == cells[0].zero_shear_viscosity_pa_s()
            });
        let profile = CarreauProfile {
            parameters,
            cells: Some(cells),
            uniform_newtonian,
        };
        let velocity =
            |i| self.accepted.q[i] / average(self.accepted.a[i], self.accepted.a[next(i, n)]);
        for i in 0..n {
            poll_at(gate, i)?;
            let strain = finite((velocity(i) - velocity(previous(i, n))) / self.spacing_m)?;
            profile.coefficients(strain, i)?;
        }
        if let Some(workspace) = &mut self.implicit_workspace {
            workspace.prepare_rheology(&profile, n, gate)?;
        }
        poll(gate)?;
        self.carreau = Some(profile);
        self.parameters.kinematic_viscosity_m2_s = reference;
        Ok(())
    }

    /// Source response at an accepted cell's actual axial velocity gradient.
    /// Returns None for a Newtonian jet; an out-of-range cell always refuses.
    /// A successful observation allocates nothing and changes neither state nor time.
    pub fn carreau_sample_at(
        &self,
        cell: usize,
    ) -> Result<Option<ThermalRheologyState>, NonlinearJetError> {
        let n = self.accepted.a.len();
        if cell >= n {
            return Err(NonlinearJetError::InvalidInput("jet rheology cell"));
        }
        let Some(parameters) = self.carreau.as_ref() else {
            return Ok(None);
        };
        let left = previous(cell, n);
        let velocity =
            |i| self.accepted.q[i] / average(self.accepted.a[i], self.accepted.a[next(i, n)]);
        let strain = finite((velocity(cell) - velocity(left)) / self.spacing_m)?;
        parameters.sample(strain, cell).map(Some)
    }
}

#[cfg(test)]
mod tests;
