//! Quasisteady viscous flow through a uniform rectangular moving slit.
//!
//! For width w, gap h, streamwise length L and mean speed v, the admitted
//! reduced law is dp = 12 mu L v / h^2 + rho v |v| / 2, with U = w h v.
//! The first term is fully developed, no-slip, parallel-plate resistance;
//! the second is the existing separated-jet loss. Explicit [`StapleRecovery`]
//! can recover part of that jet head downstream while retaining all wall loss.
//! Width is the physical slit width, not a fitted discharge coefficient.
//!
//! A P1 gap varying ACROSS the width is integrated as parallel local channels,
//! including its exact positive portion through closure. No gap is replaced by
//! area / width, which would change the h^-3 resistance. The streamwise length
//! is common; this does not resolve transverse flow, unsteady boundary layers,
//! entrance flow, pressure recovery, compressibility or turbulence.
//! Material parameters are explicit; existing inviscid callers are unchanged.
//!
//! That fully developed law is the one every constructor here builds. The reed channel of a single-reed instrument has an
//! entrance region it does not model: [`ViscousSlit::separated_jet_from_gas`] builds the separated-jet law of van Zon et al. 1990
//! instead (the module `separated_jet`, which carries the equations and the print checks), a uniform-gap law in `dp(v)` form
//! with the same bound-gas rules. `evaluate`, `characteristic_pressure` and `characteristic_junction_with_recovery` follow
//! whichever law the slit was built with; the spatial and condensed paths (`linear_gap_flow`, `fixed_opening_pressure_drop`)
//! are the fully developed law's and refuse the jet law.

mod separated_jet;
#[cfg(test)]
mod jet_tests;

use crate::acoustic_realize::AcousticRealizeError;
use fs_material::gas::GasState;
use fs_math::det;

/// Which pressure-drop law a [`ViscousSlit`] follows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SlitLaw {
    /// `dp = 12 mu L v / h^2 + rho v |v| / 2`: fully developed, no entrance region.
    FullyDeveloped,
    /// van Zon et al. 1990 (see `separated_jet`): a separated entrance that reattaches at `S = 2 h`, boundary layers, then Poiseuille.
    SeparatedJet,
}

/// Explicit fluid viscosity and uniform slit length. Both are strictly positive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViscousSlit {
    dynamic_viscosity_pa_s: f64,
    length_m: f64,
    gas: Option<GasState>,
    law: SlitLaw,
}

/// Signed flow, local tangents and nonnegative loss of the reduced slit law.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscousSlitFlow {
    /// Signed volume flow from upstream to downstream [m^3/s].
    pub flow_m3_s: f64,
    /// dU/d(dp) at fixed physical gap and width [m^3/(s Pa)].
    pub flow_per_pressure_m3_s_pa: f64,
    /// dU/dh at fixed pressure, width, viscosity and length [m^2/s].
    pub flow_per_opening_m2_s: f64,
    /// Fully developed wall-friction loss [W].
    pub viscous_dissipation_w: f64,
    /// Separated-jet loss [W], without pressure recovery.
    pub separation_dissipation_w: f64,
}

/// Integrated flow of parallel channels across a nonuniform physical slit.
/// Its single pressure drop supplies nonnegative loss `dp * flow_m3_s`.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscousProfileFlow {
    /// Signed physical volume flow [m³/s].
    pub flow_m3_s: f64,
    /// Nonnegative derivative of integrated flow with pressure [m³/(s Pa)].
    pub flow_per_pressure_m3_s_pa: f64,
}

/// Refusal text when the staple-recovery law leaves its domain. Tests match it by name.
pub const STAPLE_RECOVERY_DOMAIN: &str = "staple pressure recovery left its domain: C_P (S/S_c)^2 reached 1, so the reed opening area S is at or above S_c/sqrt(C_P) and the recovery would cancel the whole jet loss";

/// Pressure recovery in the conical staple of a double reed: Almeida, Vergez and Caussé 2007 (JASA 121(1) 536-546,
/// arXiv physics/0607011v1 pp. 8-9, eqs. 14-17). The reed opening S passes the jet; the channel then widens from the
/// constriction of area `S_c` into a cone, and a fraction `C_P` of the dynamic head at the constriction comes back as
/// pressure. With `(dp)_c = p_m - p_c` the drop across the opening and `q` the flow through it,
///
/// * flow: `q = S sqrt(2 (dp)_c / rho)` (the inviscid Bernoulli law of `moving::characteristic_junction_with_recovery`, or
///   the viscous slit law of [`ViscousSlit`] with its own linear term when the card has `slit_length`);
/// * the blade is driven by `(dp)_c` (eq. 16), not by the total mouth-to-bore drop;
/// * bore input: `p_r = p_m - (dp)_c + C_P (1/2) rho (q / S_c)^2` (eq. 17).
///
/// Both numbers are one measured reed's: `C_P = 0.8` (fitted on both pressure branches, Fig. 15) and
/// `S_c = 4e-6 m^2` (Fig. 14), a ready-to-play natural-cane oboe reed, humid air, masses added, no lips. They are borrowed
/// for any other reed. The law is a quasi-static one: `q` is the JET flow (the flow through the opening), and the swept
/// flow of the moving blade and any body flow are not carried through the staple. Recovery acts on forward flow only; a
/// reversed flow has none. It is singular where `C_P (S/S_c)^2 >= 1` (the recovery would cancel the whole jet loss and
/// the junction would create energy): there it refuses, by [`STAPLE_RECOVERY_DOMAIN`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StapleRecovery {
    c_p: f64,
    constriction_area_m2: f64,
}

impl StapleRecovery {
    /// Admit a recovery coefficient `C_P` in `[0, 1]` (0 is no recovery, 1 complete recovery) and a positive constriction area.
    ///
    /// # Errors
    /// A coefficient outside `[0, 1]` or a nonpositive or nonfinite area.
    pub fn new(c_p: f64, constriction_area_m2: f64) -> Result<Self, AcousticRealizeError> {
        if !c_p.is_finite() || !(0.0..=1.0).contains(&c_p)
            || !constriction_area_m2.is_finite() || constriction_area_m2 <= 0.0
        {
            return Err(invalid("staple recovery needs a coefficient C_P in [0, 1] and a positive finite constriction area"));
        }
        Ok(Self { c_p, constriction_area_m2 })
    }

    /// The pressure recovery coefficient `C_P` (eq. 14).
    #[must_use]
    pub const fn c_p(self) -> f64 { self.c_p }

    /// The constriction area `S_c` at the diffuser input [m^2].
    #[must_use]
    pub const fn constriction_area_m2(self) -> f64 { self.constriction_area_m2 }

    /// `1 - C_P (S/S_c)^2` for an opening area `S`: the factor on the jet's dynamic head once the recovery is in.
    ///
    /// # Errors
    /// [`STAPLE_RECOVERY_DOMAIN`] where it is not positive.
    pub fn forward_factor(self, area_m2: f64) -> Result<f64, AcousticRealizeError> {
        let ratio = area_m2 / self.constriction_area_m2;
        let factor = 1.0 - (self.c_p * ratio) * ratio;
        if !factor.is_finite() || factor <= 0.0 {
            return Err(AcousticRealizeError::Reed { what: STAPLE_RECOVERY_DOMAIN });
        }
        Ok(factor)
    }

    /// The recovered pressure `C_P (1/2) rho (q / S_c)^2` [Pa] for a forward flow, and zero for none or a reversed one.
    #[must_use]
    pub fn recovered_pressure_pa(self, density_kg_m3: f64, flow_m3_s: f64) -> f64 {
        if flow_m3_s <= 0.0 { return 0.0; }
        let velocity = flow_m3_s / self.constriction_area_m2;
        self.c_p * ((0.5 * density_kg_m3) * velocity) * velocity
    }
}

/// One same-sample characteristic junction of the slit with staple recovery.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RecoveredJunction {
    /// Outgoing characteristic wave [Pa]; the bore input pressure is this plus the incoming wave.
    pub outgoing_pa: f64,
    /// Signed flow through the reed opening [m^3/s].
    pub flow_m3_s: f64,
    /// `(dp)_c`: the drop across the opening itself (the viscous plus separated-jet law), which drives the blade [Pa].
    pub constriction_drop_pa: f64,
    /// Pressure recovered downstream of the constriction [Pa].
    pub recovered_pa: f64,
}

impl ViscousSlitFlow {
    /// Sum of the two physical losses [W]. Equal to dp * U up to rounding.
    #[must_use]
    pub fn dissipated_power_w(self) -> f64 {
        self.viscous_dissipation_w + self.separation_dissipation_w
    }
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

fn numerical() -> AcousticRealizeError {
    AcousticRealizeError::Reed {
        what: "viscous slit coefficients or result cannot be represented as finite physical values",
    }
}

impl ViscousSlit {
    /// Admit SI viscosity and streamwise slit length, without inferred defaults.
    ///
    /// # Errors
    /// Nonpositive/nonfinite inputs or an unrepresentable resistance numerator.
    pub fn new(dynamic_viscosity_pa_s: f64, length_m: f64)
        -> Result<Self, AcousticRealizeError>
    {
        if !dynamic_viscosity_pa_s.is_finite() || dynamic_viscosity_pa_s <= 0.0
            || !length_m.is_finite() || length_m <= 0.0
        {
            return Err(invalid("viscous slit requires positive finite viscosity and streamwise length"));
        }
        let numerator = 12.0 * (dynamic_viscosity_pa_s * length_m);
        if !numerator.is_finite() || numerator <= 0.0 { return Err(numerical()); }
        Ok(Self { dynamic_viscosity_pa_s, length_m, gas: None, law: SlitLaw::FullyDeveloped })
    }

    /// Bind viscosity and density to one complete material gas snapshot.
    /// Composition, pressure and temperature therefore reach the channel
    /// without an air constant; its streamwise length remains explicit geometry.
    pub fn from_gas(gas: &GasState, length_m: f64) -> Result<Self, AcousticRealizeError> {
        gas.validate_coherence().map_err(|_| invalid("viscous slit requires one coherent source gas"))?;
        let mut law = Self::new(gas.dynamic_viscosity, length_m)?;
        law.gas = Some(*gas);
        Ok(law)
    }

    /// The reed channel's separated-jet law (van Zon, Hirschberg, Gilbert and Wijnands 1990, eqs. 1-8; see the module
    /// `separated_jet`) bound to one complete material gas: viscosity and density from the gas, `length_m` the reed channel's
    /// streamwise length. Uniform-gap, stationary, laminar and 2-D, with the printed `S = 2 h`.
    ///
    /// # Errors
    /// As [`Self::from_gas`].
    pub fn separated_jet_from_gas(gas: &GasState, length_m: f64) -> Result<Self, AcousticRealizeError> {
        let mut law = Self::from_gas(gas, length_m)?;
        law.law = SlitLaw::SeparatedJet;
        Ok(law)
    }

    /// Whether this slit follows the separated-jet law (otherwise the fully developed one).
    #[must_use]
    pub const fn is_separated_jet(&self) -> bool { matches!(self.law, SlitLaw::SeparatedJet) }

    /// Original material state, when constructed from a complete source gas.
    #[must_use]
    pub const fn gas(&self) -> Option<&GasState> { self.gas.as_ref() }

    // A sourced channel follows the same new complete gas as its acoustic
    // owner. An explicitly supplied viscosity remains an explicit override.
    // The quasisteady law carries no fluid momentum or independent storage.
    pub(super) fn rebind_source_gas(self, gas: &GasState)
        -> Result<Self, AcousticRealizeError>
    {
        if self.gas.is_none() { return Ok(self); }
        let mut next = Self::from_gas(gas, self.length_m)?;
        next.law = self.law;
        Ok(next)
    }

    fn admit_density(self, density: f64) -> Result<(), AcousticRealizeError> {
        if self.gas.is_some_and(|gas| (gas.density - density).abs()
            > 64.0 * f64::EPSILON * gas.density.abs().max(density.abs())) {
            return Err(invalid("viscous slit viscosity and density must use the same source gas"));
        }
        Ok(())
    }

    /// Supplied Newtonian dynamic viscosity [Pa s].
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(self) -> f64 { self.dynamic_viscosity_pa_s }

    /// Supplied streamwise length [m].
    #[must_use]
    pub const fn length_m(self) -> f64 { self.length_m }

    fn linear_pressure(self, opening: f64) -> Result<f64, AcousticRealizeError> {
        let value = (12.0 * (self.dynamic_viscosity_pa_s * self.length_m) / opening) / opening;
        if !value.is_finite() || value <= 0.0 { return Err(numerical()); }
        Ok(value)
    }

    /// Evaluate dp = 12 mu L U / (w h^3) + rho U |U| / (2 w^2 h^2).
    /// The open-slit pressure tangent is finite even at dp = 0. A geometrically
    /// closed slit returns exactly zero flow and both zero tangents.
    ///
    /// # Errors
    /// Invalid geometry, density or pressure, or an unrepresentable result.
    pub fn evaluate(self, width_m: f64, opening_m: f64, density_kg_m3: f64, pressure_drop_pa: f64)
        -> Result<ViscousSlitFlow, AcousticRealizeError>
    {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if !pressure_drop_pa.is_finite() {
            return Err(invalid("viscous slit pressure drop must be finite"));
        }
        if opening_m == 0.0 { return Ok(ViscousSlitFlow::default()); }
        if self.is_separated_jet() { return self.jet_evaluate(width_m, area, opening_m, density_kg_m3, pressure_drop_pa); }
        let linear = self.linear_pressure(opening_m)?;
        let velocity = solve_velocity(linear, density_kg_m3, pressure_drop_pa)?;
        let denominator = linear + density_kg_m3 * velocity.abs();
        if !denominator.is_finite() || denominator <= 0.0 { return Err(numerical()); }
        let flow = area * velocity;
        // Differentiating in mean velocity avoids divisions by U at rest and
        // avoids forming R ~ h^-3 and the cancellation-prone quadratic root.
        let result = ViscousSlitFlow {
            flow_m3_s: flow,
            flow_per_pressure_m3_s_pa: area / denominator,
            flow_per_opening_m2_s: (width_m * velocity) * (1.0 + 2.0 * (linear / denominator)),
            viscous_dissipation_w: (linear * velocity) * flow,
            separation_dissipation_w: ((0.5 * density_kg_m3 * velocity) * velocity.abs()) * flow,
        };
        if ![result.flow_m3_s, result.flow_per_pressure_m3_s_pa,
            result.flow_per_opening_m2_s, result.viscous_dissipation_w,
            result.separation_dissipation_w, result.dissipated_power_w()]
            .iter().all(|v| v.is_finite())
        {
            return Err(numerical());
        }
        Ok(result)
    }

    /// Solve the SAME-sample characteristic junction, including supplied swept
    /// and body flow, rather than filtering an inviscid result afterward.
    ///
    /// With A = w h and d = Pmouth - 2 Pminus - Z Umoving, solve
    /// (12 mu L / h^2 + Z A) v + rho v |v| / 2 = d. Return Pplus [Pa].
    /// No nonlinear iteration, state mutation or allocation occurs here.
    ///
    /// # Errors
    /// Invalid physical inputs, nonfinite intermediate pressure or overflow.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_pressure(
        self, width_m: f64, opening_m: f64, density_kg_m3: f64,
        impedance_pa_s_m3: f64, incoming_pressure_pa: f64,
        upstream_pressure_pa: f64, moving_flow_m3_s: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if ![impedance_pa_s_m3, incoming_pressure_pa, upstream_pressure_pa, moving_flow_m3_s]
            .iter().all(|v| v.is_finite()) || impedance_pa_s_m3 <= 0.0
        {
            return Err(invalid("viscous slit junction requires finite pressures/flow and positive impedance"));
        }
        let no_jet = incoming_pressure_pa + impedance_pa_s_m3 * moving_flow_m3_s;
        let no_drop = upstream_pressure_pa - incoming_pressure_pa;
        let drive = no_drop - no_jet;
        if ![no_jet, no_drop, drive].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        if opening_m == 0.0 || drive == 0.0 { return Ok(no_jet); }
        if self.is_separated_jet() {
            let jet = self.jet_velocity(opening_m, density_kg_m3, impedance_pa_s_m3 * area, 0.0, drive)?;
            return jet_outgoing(no_jet, no_drop, impedance_pa_s_m3 * area, jet);
        }
        let linear = self.linear_pressure(opening_m)?;
        let wave_load = impedance_pa_s_m3 * area;
        let velocity = solve_velocity(linear + wave_load, density_kg_m3, drive)?;
        let jet_pressure = wave_load * velocity;
        // Keep the constitutive drop: reconstructing it as drive - Z U would
        // erase small pressure drops when the wave load dominates resistance.
        let dp = linear * velocity + (0.5 * density_kg_m3 * velocity) * velocity.abs();
        let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
            no_drop - dp
        } else {
            no_jet + jet_pressure
        };
        if ![dp, jet_pressure, outgoing].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        Ok(outgoing)
    }

    /// [`Self::characteristic_pressure`] with the staple's pressure recovery ([`StapleRecovery`]): the SAME-sample
    /// junction `(dp)_c(q) - C_P (1/2) rho (q/S_c)^2 + Z q = d` for forward flow, with `d = Pmouth - 2 Pminus - Z Umoving`
    /// and `(dp)_c(q) = linear U/A + rho (U/A)|U/A| / 2` the slit law. For forward flow that is the slit's own quadratic with
    /// the density scaled by `1 - C_P (S/S_c)^2`; reversed flow has no recovery and is exactly the original junction.
    /// Returns the outgoing wave and the pieces the blade and the energy account need.
    ///
    /// # Errors
    /// As [`Self::characteristic_pressure`], and [`STAPLE_RECOVERY_DOMAIN`] where forward flow meets `C_P (S/S_c)^2 >= 1`.
    #[allow(clippy::too_many_arguments)]
    pub fn characteristic_junction_with_recovery(
        self, recovery: StapleRecovery, width_m: f64, opening_m: f64, density_kg_m3: f64,
        impedance_pa_s_m3: f64, incoming_pressure_pa: f64,
        upstream_pressure_pa: f64, moving_flow_m3_s: f64,
    ) -> Result<RecoveredJunction, AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if ![impedance_pa_s_m3, incoming_pressure_pa, upstream_pressure_pa, moving_flow_m3_s]
            .iter().all(|v| v.is_finite()) || impedance_pa_s_m3 <= 0.0
        {
            return Err(invalid("viscous slit junction requires finite pressures/flow and positive impedance"));
        }
        let no_jet = incoming_pressure_pa + impedance_pa_s_m3 * moving_flow_m3_s;
        let no_drop = upstream_pressure_pa - incoming_pressure_pa;
        let drive = no_drop - no_jet;
        if ![no_jet, no_drop, drive].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        if opening_m == 0.0 || drive == 0.0 {
            // No flow, so no recovery: the whole mouth-to-bore difference `drive` falls across the shut reed, and it is
            // that difference the blade feels (zero only when there is no difference).
            return Ok(RecoveredJunction {
                outgoing_pa: no_jet, flow_m3_s: 0.0, constriction_drop_pa: drive, recovered_pa: 0.0,
            });
        }
        if self.is_separated_jet() {
            return self.jet_junction_with_recovery(recovery, area, opening_m, density_kg_m3, impedance_pa_s_m3,
                no_jet, no_drop, drive);
        }
        let linear = self.linear_pressure(opening_m)?;
        let wave_load = impedance_pa_s_m3 * area;
        // The junction is monotone in the mean velocity and zero at rest, so forward flow is a positive drive.
        let jet_density = if drive > 0.0 {
            density_kg_m3 * recovery.forward_factor(area)?
        } else {
            density_kg_m3
        };
        let velocity = solve_velocity(linear + wave_load, jet_density, drive)?;
        let flow = area * velocity;
        let jet_pressure = wave_load * velocity;
        // The constitutive drops, kept as such: reconstructing them from the drive would erase small drops.
        let constriction_drop = linear * velocity + (0.5 * density_kg_m3 * velocity) * velocity.abs();
        let recovered = recovery.recovered_pressure_pa(density_kg_m3, flow);
        let dp = constriction_drop - recovered;
        let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
            no_drop - dp
        } else {
            no_jet + jet_pressure
        };
        if ![flow, constriction_drop, recovered, dp, outgoing].iter().all(|v| v.is_finite()) {
            return Err(numerical());
        }
        Ok(RecoveredJunction {
            outgoing_pa: outgoing, flow_m3_s: flow, constriction_drop_pa: constriction_drop,
            recovered_pa: recovered,
        })
    }

    /// Eliminate dp + Z U(dp) = drive after condensing mechanical swept flow.
    /// Return dp directly: subtracting large wave pressures can erase it.
    pub(super) fn fixed_opening_pressure_drop(
        self, width_m: f64, opening_m: f64, density_kg_m3: f64,
        impedance_pa_s_m3: f64, drive_pa: f64,
    ) -> Result<f64, AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 || !drive_pa.is_finite() {
            return Err(invalid("viscous slit condensation requires positive impedance and finite drive"));
        }
        if opening_m == 0.0 || drive_pa == 0.0 { return Ok(drive_pa); }
        if self.is_separated_jet() { return Err(jet_refusal()); }
        let linear = self.linear_pressure(opening_m)?;
        let velocity = solve_velocity(linear + impedance_pa_s_m3 * area, density_kg_m3, drive_pa)?;
        let dp = linear * velocity + (0.5 * density_kg_m3 * velocity) * velocity.abs();
        if !dp.is_finite() { return Err(numerical()); }
        Ok(dp)
    }

    /// The same fully developed slit after mechanical swept-flow condensation,
    /// including downstream staple recovery. Returns (blade drop, total drop).
    ///
    /// With D = 1 + Z * pressure_mobility, the equation is
    /// `p_blade(v) - recovery(v)/D + (Z/D) A v = drive`.
    /// Only the quadratic jet head is recovered; the linear no-slip wall loss
    /// remains intact. This is a quasisteady uniform slit, not a spatial or
    /// separated-entrance model. The original passive recovery domain applies.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn fixed_opening_pressure_drop_with_recovery(
        self, recovery: StapleRecovery, width_m: f64, opening_m: f64,
        density_kg_m3: f64, impedance_pa_s_m3: f64, denominator: f64, drive_pa: f64,
    ) -> Result<(f64, f64), AcousticRealizeError> {
        let area = admitted_area(width_m, opening_m, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0
            || !denominator.is_finite() || denominator < 1.0 || !drive_pa.is_finite() {
            return Err(invalid("viscous slit recovery condensation requires positive impedance, finite drive and denominator at least one"));
        }
        if self.is_separated_jet() { return Err(jet_refusal()); }
        // Retain original arithmetic when no pressure can be recovered.
        if recovery.c_p() == 0.0 || drive_pa <= 0.0 || opening_m == 0.0 {
            let dp = self.fixed_opening_pressure_drop(
                width_m, opening_m, density_kg_m3, impedance_pa_s_m3, drive_pa)?;
            return Ok((dp, dp));
        }
        let remaining_head = recovery.forward_factor(area)?;
        let condensed_head = 1.0 - (1.0 - remaining_head) / denominator;
        let linear = self.linear_pressure(opening_m)?;
        let velocity = solve_velocity(
            linear + impedance_pa_s_m3 * area, density_kg_m3 * condensed_head, drive_pa)?;
        let wall_drop = linear * velocity;
        let jet_drop = (0.5 * density_kg_m3 * velocity) * velocity;
        let blade = wall_drop + jet_drop;
        // Do not recover viscous wall friction or subtract two nearly equal
        // full drops when the allowed recovery approaches the jet head.
        let total = wall_drop + remaining_head * jet_drop;
        if ![blade, total].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        Ok((blade, total))
    }

    /// Integrate this same law over a linearly varying widthwise gap. Negative
    /// gaps contribute no flow; a partially closed edge retains its actual open
    /// fraction. The pressure derivative is integrated analytically as well.
    ///
    /// With K=12 mu L and C=2 rho |dp|, flow per width is
    /// `(sqrt(K²+C h⁴)-K)/(rho h)`. Its primitive on h>0 is
    /// `[r-K ln(K+r)]/(2 rho)`, r=sqrt(K²+C h⁴). The divided difference below
    /// rationalizes delta-r and log1p so nearly equal gaps and the laminar limit
    /// do not subtract nearly equal primitives. This is exact for the declared
    /// parallel-channel P1 model, not a spatial Navier--Stokes solution.
    pub fn linear_gap_flow(self, width_m: f64, gaps_m: [f64; 2], density_kg_m3: f64,
        pressure_drop_pa: f64) -> Result<ViscousProfileFlow, AcousticRealizeError>
    {
        if gaps_m.iter().any(|g| !g.is_finite()) || !pressure_drop_pa.is_finite() {
            return Err(invalid("viscous slit profile requires finite gaps and pressure"));
        }
        let lower = gaps_m[0].min(gaps_m[1]);
        let upper = gaps_m[0].max(gaps_m[1]).max(0.0);
        admitted_area(width_m, upper, density_kg_m3)?;
        self.admit_density(density_kg_m3)?;
        if upper == 0.0 { return Ok(ViscousProfileFlow::default()); }
        if self.is_separated_jet() { return Err(jet_refusal()); }
        if lower == upper {
            let flow = self.evaluate(width_m, upper, density_kg_m3, pressure_drop_pa)?;
            return Ok(ViscousProfileFlow { flow_m3_s: flow.flow_m3_s,
                flow_per_pressure_m3_s_pa: flow.flow_per_pressure_m3_s_pa });
        }
        let open_fraction = if lower >= 0.0 { 1.0 } else {
            let scale = upper.max(-lower);
            (upper / scale) / (upper / scale + (-lower) / scale)
        };
        let ratio = lower.max(0.0) / upper;
        let linear = self.linear_pressure(upper)?;
        let quadratic = (std::f64::consts::SQRT_2 * det::sqrt(density_kg_m3))
            * det::sqrt(pressure_drop_pa.abs());
        if !quadratic.is_finite() { return Err(numerical()); }
        let scale = linear.max(quadratic);
        let (l, s) = (linear / scale, quadratic / scale);
        let r0 = det::hypot(l, s * ratio * ratio);
        let r1 = det::hypot(l, s);
        let polynomial = (1.0 + ratio) * (1.0 + ratio * ratio);
        let tangent = (width_m * open_fraction) * (upper / scale)
            * (polynomial / (2.0 * (r0 + r1)));
        let delta = s * s * (1.0 - ratio) * polynomial / (r0 + r1);
        let log_argument = delta / (l + r0);
        let log_ratio = if !log_argument.is_finite() { 0.0 } else {
            let sum = 1.0 + log_argument;
            if sum == 1.0 { 1.0 } else { det::ln(sum) / (sum - 1.0) }
        };
        let factor = if l == 0.0 { 1.0 } else { 1.0 - (l / (l + r0)) * log_ratio };
        let flow = (pressure_drop_pa * tangent) * (2.0 * factor);
        if !flow.is_finite() || !tangent.is_finite() || tangent < 0.0 {
            return Err(numerical());
        }
        Ok(ViscousProfileFlow { flow_m3_s: flow, flow_per_pressure_m3_s_pa: tangent })
    }
}

fn jet_refusal() -> AcousticRealizeError {
    invalid("the separated-jet slit law is defined on a uniform reed gap only (spatial and condensed paths use the fully developed law)")
}

/// The speed solved for a jet-law junction, with the slit's own signed drop at it.
#[derive(Debug, Clone, Copy)]
struct JetSolution {
    /// Signed mean speed [m/s].
    velocity: f64,
    /// Signed constitutive drop of the slit alone [Pa].
    drop: f64,
    /// The law at that speed (unsigned speed), with its tangents.
    law: separated_jet::Drop,
}

/// The outgoing characteristic wave from the same two candidates the fully developed junction picks between.
fn jet_outgoing(no_jet: f64, no_drop: f64, wave_load: f64, jet: JetSolution) -> Result<f64, AcousticRealizeError> {
    let jet_pressure = wave_load * jet.velocity;
    let outgoing = if no_drop.abs().max(jet.drop.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
        no_drop - jet.drop
    } else {
        no_jet + jet_pressure
    };
    if ![jet.drop, jet_pressure, outgoing].iter().all(|v| v.is_finite()) { return Err(numerical()); }
    Ok(outgoing)
}

impl ViscousSlit {
    /// Signed speed and slit drop of `dp(v) + wave_load v - recovery v^2 = drive` (the recovery acts on forward flow only, so the
    /// caller passes zero for a reversed drive).
    fn jet_velocity(self, opening_m: f64, density_kg_m3: f64, wave_load: f64, recovery: f64, drive: f64)
        -> Result<JetSolution, AcousticRealizeError>
    {
        let speed = separated_jet::speed_for_drive(
            opening_m, self.length_m, density_kg_m3, self.dynamic_viscosity_pa_s, wave_load, recovery, drive.abs(),
        ).ok_or_else(numerical)?;
        let sign = if drive < 0.0 { -1.0 } else { 1.0 };
        let law = separated_jet::drop_at(speed, opening_m, self.length_m, density_kg_m3, self.dynamic_viscosity_pa_s);
        if !law.dp.is_finite() { return Err(numerical()); }
        Ok(JetSolution { velocity: sign * speed, drop: sign * law.dp, law })
    }

    /// [`Self::evaluate`] on the separated-jet law. At rest on the contracted branch (`L / h < 3`, `dp` quadratic in `v`) the
    /// pressure tangent is infinite and reported as such; every other field is finite.
    fn jet_evaluate(self, width_m: f64, area: f64, opening_m: f64, density_kg_m3: f64, pressure_drop_pa: f64)
        -> Result<ViscousSlitFlow, AcousticRealizeError>
    {
        let jet = self.jet_velocity(opening_m, density_kg_m3, 0.0, 0.0, pressure_drop_pa)?;
        let speed = jet.velocity.abs();
        let sign = if jet.velocity < 0.0 { -1.0 } else { 1.0 };
        let law = jet.law;
        let flow = area * jet.velocity;
        let tangent = if law.d_dv > 0.0 { area / law.d_dv } else { f64::INFINITY };
        // U = w h v(h, dp): dU/dh = w v - sign w h (d dp/dh) / (d dp/dv) (both zero at rest).
        let opening_term = if law.d_dv > 0.0 { sign * (width_m * opening_m) * (law.d_dh / law.d_dv) } else { 0.0 };
        let separation = (((0.5 * density_kg_m3) * speed) * speed) * (area * speed);
        let total = law.dp * (area * speed);
        let result = ViscousSlitFlow {
            flow_m3_s: flow,
            flow_per_pressure_m3_s_pa: tangent,
            flow_per_opening_m2_s: width_m * jet.velocity - opening_term,
            // Everything beyond the mean-flow dynamic head: the entrance region, the contraction below L / h = 3, wall friction.
            viscous_dissipation_w: (total - separation).max(0.0),
            separation_dissipation_w: separation,
        };
        if ![result.flow_m3_s, result.flow_per_opening_m2_s, result.viscous_dissipation_w, result.separation_dissipation_w]
            .iter().all(|v| v.is_finite()) || result.flow_per_pressure_m3_s_pa.is_nan()
        {
            return Err(numerical());
        }
        Ok(result)
    }

    /// [`Self::characteristic_junction_with_recovery`] on the separated-jet law: the staple's recovery is `recovery * v^2` with
    /// `recovery = C_P rho (A / S_c)^2 / 2`, `rho / 2` less than the jet head it reduces (`forward_factor` is that margin).
    #[allow(clippy::too_many_arguments)]
    fn jet_junction_with_recovery(
        self, recovery: StapleRecovery, area: f64, opening_m: f64, density_kg_m3: f64, impedance_pa_s_m3: f64,
        no_jet: f64, no_drop: f64, drive: f64,
    ) -> Result<RecoveredJunction, AcousticRealizeError> {
        let wave_load = impedance_pa_s_m3 * area;
        let head_removed = if drive > 0.0 {
            0.5 * density_kg_m3 * (1.0 - recovery.forward_factor(area)?)
        } else {
            0.0
        };
        let jet = self.jet_velocity(opening_m, density_kg_m3, wave_load, head_removed, drive)?;
        let flow = area * jet.velocity;
        let recovered = recovery.recovered_pressure_pa(density_kg_m3, flow);
        let dp = jet.drop - recovered;
        let jet_pressure = wave_load * jet.velocity;
        let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
            no_drop - dp
        } else {
            no_jet + jet_pressure
        };
        if ![flow, jet.drop, recovered, dp, outgoing].iter().all(|v| v.is_finite()) { return Err(numerical()); }
        Ok(RecoveredJunction {
            outgoing_pa: outgoing, flow_m3_s: flow, constriction_drop_pa: jet.drop, recovered_pa: recovered,
        })
    }
}

/// Shared monotone passive pressure solve for a supplied spatial slit integral.
/// The constitutive law is odd and increasing. All 96 bounded evaluations are
/// at the same mechanical trial; no accepted state or pressure history changes.
pub(super) fn profiled_pressure_drop(
    impedance: f64, drive: f64,
    flow: impl Fn(f64) -> Result<ViscousProfileFlow, AcousticRealizeError>,
) -> Result<f64, AcousticRealizeError> {
    if !impedance.is_finite() || impedance <= 0.0 || !drive.is_finite() {
        return Err(invalid("profiled slit requires finite drive and positive passive impedance"));
    }
    if drive == 0.0 { return Ok(0.0); }
    let source = drive.abs();
    let (mut lo, mut hi, mut p) = (0.0, source, source);
    for _ in 0..96 {
        let value = flow(p)?;
        let load = impedance * value.flow_m3_s;
        let derivative = 1.0 + impedance * value.flow_per_pressure_m3_s_pa;
        let residual = p + load - source;
        if ![load, derivative, residual].iter().all(|x| x.is_finite())
            || value.flow_m3_s < 0.0 || value.flow_per_pressure_m3_s_pa < 0.0 {
            return Err(numerical());
        }
        if residual.abs() <= 32.0 * f64::EPSILON * source {
            return Ok(drive.signum() * p);
        }
        if residual > 0.0 { hi = p; } else { lo = p; }
        let newton = p - residual / derivative;
        p = if newton > lo && newton < hi { newton } else { f64::midpoint(lo, hi) };
        if p == lo || p == hi { break; }
    }
    Err(invalid("profiled viscous slit pressure balance did not converge within its bounded work"))
}

/// Same spatial pressure solve before mechanical swept-flow condensation.
pub(super) fn profiled_characteristic_pressure(
    impedance: f64, incoming: f64, upstream: f64, moving: f64,
    flow: impl Fn(f64) -> Result<ViscousProfileFlow, AcousticRealizeError>,
) -> Result<f64, AcousticRealizeError> {
    let no_jet = incoming + impedance * moving;
    let no_drop = upstream - incoming;
    let dp = profiled_pressure_drop(impedance, no_drop - no_jet, &flow)?;
    let jet_pressure = impedance * flow(dp)?.flow_m3_s;
    let outgoing = if no_drop.abs().max(dp.abs()) <= no_jet.abs().max(jet_pressure.abs()) {
        no_drop - dp
    } else { no_jet + jet_pressure };
    if !outgoing.is_finite() { return Err(numerical()); }
    Ok(outgoing)
}

fn admitted_area(width: f64, opening: f64, rho: f64) -> Result<f64, AcousticRealizeError> {
    if ![width, opening, rho].iter().all(|v| v.is_finite())
        || width <= 0.0 || opening < 0.0 || rho <= 0.0
    {
        return Err(invalid("viscous slit requires positive finite width/density and a nonnegative finite gap"));
    }
    let area = width * opening;
    if !area.is_finite() || (opening > 0.0 && area == 0.0) { return Err(numerical()); }
    Ok(area)
}

fn solve_velocity(linear: f64, rho: f64, drive: f64) -> Result<f64, AcousticRealizeError> {
    if !linear.is_finite() || linear <= 0.0 { return Err(numerical()); }
    if drive == 0.0 { return Ok(0.0); }
    // v = 2 |d| / (a + sqrt(a^2 + 2 rho |d|)), rationalized and scaled.
    // sqrt(|d|) keeps sub-threshold pressures continuous; a is never squared.
    let quadratic = (std::f64::consts::SQRT_2 * det::sqrt(rho)) * det::sqrt(drive.abs());
    if !quadratic.is_finite() { return Err(numerical()); }
    let scale = linear.max(quadratic);
    let an = linear / scale;
    let bn = quadratic / scale;
    let denominator = 0.5 * (an + det::sqrt(an * an + bn * bn));
    let velocity = drive.signum() * ((drive.abs() / scale) / denominator);
    if !velocity.is_finite() { return Err(numerical()); }
    Ok(velocity)
}

#[cfg(test)]
#[allow(clippy::float_cmp)] // Exact closure, rest, symmetry and refusal pins.
mod tests {
    use super::*;

    fn law() -> ViscousSlit { ViscousSlit::new(1.81e-5, 1e-3).unwrap() }

    fn close(actual: f64, expected: f64, relative: f64) {
        assert!((actual - expected).abs() <= relative * expected.abs().max(1e-250),
            "actual={actual:e}, expected={expected:e}");
    }

    #[test]
    fn g3_material_gas_supplies_viscosity_and_cannot_be_mixed_with_another_density() {
        use fs_material::gas::GasSpec;
        let cold = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let hot = GasState::try_new(&GasSpec::dry_air_ussa1976(), 393.15, 101325.0).unwrap();
        let a = ViscousSlit::from_gas(&cold, 0.001).unwrap();
        let b = ViscousSlit::from_gas(&hot, 0.001).unwrap();
        assert_eq!(a.gas(), Some(&cold));
        assert_eq!(a.dynamic_viscosity_pa_s(), cold.dynamic_viscosity);
        let flow = |law: ViscousSlit, gas: GasState| law.evaluate(0.01, 1e-5, gas.density, 100.0).unwrap().flow_m3_s;
        assert!(flow(b, hot) < flow(a, cold), "a thin channel must respond to material viscosity, not only density");
        assert!(a.evaluate(0.01, 1e-5, hot.density, 100.0).is_err());
        let mut corrupt = cold; corrupt.dynamic_viscosity *= 2.0;
        assert!(ViscousSlit::from_gas(&corrupt, 0.001).is_err());
    }

    #[test]
    fn g1_spatial_laminar_limit_integrates_cubic_gap_and_partial_closure() {
        let (width, h, dp) = (0.013, 4e-5, 1e-8);
        let k = 12.0 * 1.81e-5 * 1e-3;
        let triangular = law().linear_gap_flow(width, [0.0, 2.0*h], 1.2, dp).unwrap();
        let mean_gap = law().evaluate(width, h, 1.2, dp).unwrap();
        close(triangular.flow_m3_s, width*dp*(2.0*h).powi(3)/(4.0*k), 1e-10);
        close(triangular.flow_m3_s/mean_gap.flow_m3_s, 2.0, 1e-10);
        let partial = law().linear_gap_flow(width, [-2.0*h, 2.0*h], 1.2, dp).unwrap();
        close(partial.flow_m3_s, 0.5*triangular.flow_m3_s, 1e-14);
        let tangent = law().linear_gap_flow(width, [0.0, 2.0*h], 1.2, 0.0).unwrap();
        close(tangent.flow_per_pressure_m3_s_pa, width*(2.0*h).powi(3)/(4.0*k), 1e-14);
        assert_eq!(law().linear_gap_flow(width, [-h, 0.0], 1.2, dp).unwrap(), ViscousProfileFlow::default());
    }

    #[test]
    fn g1_spatial_integral_matches_independent_point_quadrature_and_pressure_tangent() {
        let (width, rho) = (0.013, 1.2);
        for gaps in [[-4e-5, 8e-5], [1e-5, 1e-4], [1e-4, 1e-4*(1.0+1e-14)]] {
            for dp in [-2000.0, 0.1, 2000.0] {
                let exact = law().linear_gap_flow(width, gaps, rho, dp).unwrap();
                let mut numerical = 0.0;
                for i in 0..8192 {
                    let t = (i as f64 + 0.5)/8192.0;
                    let h = ((1.0-t)*gaps[0]+t*gaps[1]).max(0.0);
                    numerical += law().evaluate(width/8192.0, h, rho, dp).unwrap().flow_m3_s;
                }
                close(exact.flow_m3_s, numerical, 4e-8);
                let delta = dp.abs()*1e-5;
                let left = law().linear_gap_flow(width, gaps, rho, dp-delta).unwrap();
                let right = law().linear_gap_flow(width, gaps, rho, dp+delta).unwrap();
                close(exact.flow_per_pressure_m3_s_pa, (right.flow_m3_s-left.flow_m3_s)/(2.0*delta), 2e-9);
                let mid = f64::midpoint(gaps[0], gaps[1]);
                let a = law().linear_gap_flow(width/2.0, [gaps[0],mid], rho, dp).unwrap();
                let b = law().linear_gap_flow(width/2.0, [mid,gaps[1]], rho, dp).unwrap();
                close(exact.flow_m3_s, a.flow_m3_s+b.flow_m3_s, 2e-13);
                let z = 1e6;
                let solved = profiled_pressure_drop(z, dp+z*exact.flow_m3_s,
                    |p| law().linear_gap_flow(width, gaps, rho, p)).unwrap();
                close(solved, dp, 1e-12);
            }
        }
    }

    #[test]
    fn g1_manufactured_signed_flows_recover_both_pressure_losses() {
        let (width, rho) = (0.013, 1.2);
        for h in [1e-7, 1e-5, 4e-4] {
            for velocity in [-20.0_f64, -0.25, 0.0, 0.25, 20.0] {
                let pv = 12.0 * 1.81e-5 * 1e-3 * velocity / (h * h);
                let ps = 0.5 * rho * velocity * velocity.abs();
                let actual = law().evaluate(width, h, rho, pv + ps).unwrap();
                let expected = width * h * velocity;
                close(actual.flow_m3_s, expected, 3e-14);
                close(actual.viscous_dissipation_w, pv * expected, 6e-14);
                close(actual.separation_dissipation_w, ps * expected, 6e-14);
                close(actual.dissipated_power_w(), (pv + ps) * expected, 6e-14);
                assert!(actual.viscous_dissipation_w >= 0.0);
                assert!(actual.separation_dissipation_w >= 0.0);
            }
        }
    }

    #[test]
    fn g1_laminar_limit_has_finite_rest_tangent_and_cubic_gap_scaling() {
        let (width, h, rho) = (0.013, 1e-5, 1.2);
        let at_rest = law().evaluate(width, h, rho, 0.0).unwrap();
        close(at_rest.flow_per_pressure_m3_s_pa, width * h * h * h / (12.0 * 1.81e-5 * 1e-3), 3e-14);
        assert_eq!(at_rest.flow_m3_s, 0.0);
        assert_eq!(at_rest.flow_per_opening_m2_s, 0.0);
        let first = law().evaluate(width, h, rho, 1e-6).unwrap();
        let doubled = law().evaluate(width, 2.0 * h, rho, 1e-6).unwrap();
        close(doubled.flow_m3_s / first.flow_m3_s, 8.0, 1e-8);
    }

    #[test]
    fn g1_low_viscosity_recovers_the_separated_jet_limit() {
        let (width, h, rho, dp) = (0.013, 4e-4, 1.2_f64, 250.0_f64);
        let inviscid = width * h * (2.0 * dp / rho).sqrt();
        let nearly_inviscid = ViscousSlit::new(1e-12, 1e-5).unwrap().evaluate(width, h, rho, dp).unwrap();
        close(nearly_inviscid.flow_m3_s, inviscid, 1e-8);
        let viscous = law().evaluate(width, h, rho, dp).unwrap();
        assert!(viscous.flow_m3_s > 0.0 && viscous.flow_m3_s < inviscid);
    }

    #[test]
    fn g3_pressure_and_gap_tangents_match_independent_differences() {
        let (width, h, rho) = (0.013, 4e-5, 1.2);
        for dp in [-250.0_f64, -1e-5, 1e-5, 250.0] {
            let value = law().evaluate(width, h, rho, dp).unwrap();
            let ep = dp.abs() * 1e-5;
            let eh = h * 1e-5;
            let flow = |gap, pressure| law().evaluate(width, gap, rho, pressure).unwrap().flow_m3_s;
            close(value.flow_per_pressure_m3_s_pa, (flow(h, dp + ep) - flow(h, dp - ep)) / (2.0 * ep), 2e-8);
            close(value.flow_per_opening_m2_s, (flow(h + eh, dp) - flow(h - eh, dp)) / (2.0 * eh), 2e-8);
        }
    }

    #[test]
    fn g1_characteristic_junction_recovers_manufactured_bidirectional_flow() {
        let (width, rho, impedance) = (0.013, 1.2, 1e6);
        for h in [1e-7, 4e-5, 4e-4] {
            for velocity in [-3.0_f64, -1e-5, 1e-5, 3.0] {
                let dp = 12.0 * 1.81e-5 * 1e-3 * velocity / (h * h)
                    + 0.5 * rho * velocity * velocity.abs();
                for incoming in [-200.0, 0.0, 200.0] {
                    for moving in [-1e-6, 0.0, 1e-6] {
                        let expected = incoming + impedance * (width * h * velocity + moving);
                        let upstream = dp + expected + incoming;
                        let actual = law().characteristic_pressure(width, h, rho, impedance, incoming, upstream, moving).unwrap();
                        assert!((actual - expected).abs() <= 3e-12 * (1.0 + expected.abs()));
                    }
                }
            }
        }
    }

    #[test]
    fn g0_exact_closure_and_arbitrarily_small_signed_drives() {
        let closed = law().evaluate(0.013, 0.0, 1.2, 500.0).unwrap();
        assert_eq!(closed, ViscousSlitFlow::default());
        assert_eq!(law().characteristic_pressure(0.013, 0.0, 1.2, 1e6, 75.0, 1000.0, 2e-7).unwrap(), 75.2);
        for dp in [1e-20, 1e-100, 1e-200] {
            let forward = law().evaluate(0.013, 1e-4, 1.2, dp).unwrap();
            let reverse = law().evaluate(0.013, 1e-4, 1.2, -dp).unwrap();
            assert!(forward.flow_m3_s > 0.0);
            assert_eq!(reverse.flow_m3_s, -forward.flow_m3_s);
            assert_eq!(reverse.flow_per_pressure_m3_s_pa, forward.flow_per_pressure_m3_s_pa);
        }
    }

    #[test]
    fn g0_invalid_and_unrepresentable_inputs_are_refused() {
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(ViscousSlit::new(bad, 1e-3).is_err());
            assert!(ViscousSlit::new(1.81e-5, bad).is_err());
            assert!(law().evaluate(bad, 1e-4, 1.2, 100.0).is_err());
            assert!(law().evaluate(0.013, 1e-4, bad, 100.0).is_err());
            assert!(law().characteristic_pressure(0.013, 1e-4, 1.2, bad, 0.0, 100.0, 0.0).is_err());
        }
        for bad in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(law().evaluate(0.013, bad, 1.2, 100.0).is_err());
        }
        assert!(law().evaluate(0.013, 0.0, 1.2, f64::NAN).is_err());
        assert!(law().evaluate(0.013, 1e-200, 1.2, 100.0).is_err());
        assert!(ViscousSlit::new(f64::MAX, f64::MAX).is_err());
    }

    // ---- staple pressure recovery (Almeida, Vergez and Caussé 2007, eqs. 14-17) -----------------------------------
    // Fig. 13, increasing branch: p_M = 35 kPa, k_s = 1.04e10 kg m^-3 s^-2. Fig. 14 and Fig. 15: S_c = 4e-6 m^2, C_P = 0.8.
    const P_M: f64 = 35e3;
    const K_S: f64 = 1.04e10;
    const S_C: f64 = 4e-6;
    const C_P: f64 = 0.8;

    /// The junction with the viscous term and the bore impedance negligible, so it is the printed inviscid curve: the
    /// bore input pressure is then zero and the total drop is the mouth pressure. Returns (total drop, flow, constriction drop).
    fn quasistatic(recovery: StapleRecovery, area: f64, mouth: f64) -> (f64, f64, f64) {
        let slit = ViscousSlit::new(1e-12, 1e-9).unwrap();
        let width = 0.0056;
        let j = slit.characteristic_junction_with_recovery(
            recovery, width, area / width, 1.2, 1.0, 0.0, mouth, 0.0).unwrap();
        (mouth - j.outgoing_pa, j.flow_m3_s, j.constriction_drop_pa)
    }

    /// One point of the printed characteristic: the constriction drop `x` sets the opening `S = S_0 - x/k_s` (eq. 16); the
    /// junction is solved for the total drop at which its own constriction drop is `x`. Returns (total drop, flow).
    fn characteristic_point(recovery: StapleRecovery, x: f64) -> (f64, f64) {
        let area = P_M / K_S - x / K_S;
        let (mut lo, mut hi) = (0.0, x);
        for _ in 0..90 {
            let mid = 0.5 * (lo + hi);
            if quasistatic(recovery, area, mid).2 < x { lo = mid; } else { hi = mid; }
        }
        let (total, flow, _) = quasistatic(recovery, area, 0.5 * (lo + hi));
        (total, flow)
    }

    /// Golden-section maximum of the flow over the constriction drop. Returns (x, flow).
    fn flow_maximum(recovery: StapleRecovery) -> (f64, f64) {
        let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
        let f = |x: f64| characteristic_point(recovery, x).1;
        let (mut a, mut b) = (0.15 * P_M, 0.6 * P_M);
        let (mut c, mut d) = (b - phi * (b - a), a + phi * (b - a));
        let (mut fc, mut fd) = (f(c), f(d));
        for _ in 0..80 {
            if fc > fd {
                b = d; d = c; fd = fc;
                c = b - phi * (b - a); fc = f(c);
            } else {
                a = c; c = d; fc = fd;
                d = a + phi * (b - a); fd = f(d);
            }
        }
        let x = 0.5 * (a + b);
        (x, f(x))
    }

    #[test]
    fn g1_staple_recovery_puts_the_flow_maximum_where_the_printed_model_says() {
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        let s0 = P_M / K_S; // p_M = k_s S_0
        // Closed form from eqs. 15-17: q^2 = S^2 2x/rho with S = S_0 - x/k_s is largest at x = p_M/3, S = 2 S_0/3, and the
        // total drop there is x (1 - C_P S^2/S_c^2), because (1/2) rho (q/S_c)^2 = x S^2/S_c^2. Nothing here is a typed ratio.
        let x_star = P_M / 3.0;
        let s_star = 2.0 * s0 / 3.0;
        let expected_drop = x_star * (1.0 - C_P * (s_star / S_C) * (s_star / S_C));
        let expected_flow = s_star * (2.0 * x_star / 1.2).sqrt();
        let (x, flow) = flow_maximum(recovery);
        let (drop, _) = characteristic_point(recovery, x);
        // The recovery moves where the maximum sits on the total-drop axis, not how large it is.
        close(flow, expected_flow, 1e-6);
        close(drop, expected_drop, 1e-5);
        assert!(drop < 0.9 * x_star, "the recovery must pull the maximum below the elementary model's");
        // Without recovery the maximum is the elementary model's, at p_M/3.
        let (x0, flow0) = flow_maximum(StapleRecovery::new(0.0, S_C).unwrap());
        close(x0, x_star, 1e-5);
        close(flow0, expected_flow, 1e-6);
    }

    #[test]
    fn g0_staple_recovery_refuses_where_the_jet_loss_would_be_cancelled() {
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        let width = 0.0056;
        let limit = S_C / C_P.sqrt(); // C_P (S/S_c)^2 = 1
        let junction = |area: f64, mouth: f64| law().characteristic_junction_with_recovery(
            recovery, width, area / width, 1.2, 1e7, 0.0, mouth, 0.0);
        assert!(junction(limit * (1.0 - 1e-9), 500.0).is_ok());
        for area in [limit * (1.0 + 1e-9), 1.5 * limit] {
            match junction(area, 500.0) {
                Err(AcousticRealizeError::Reed { what }) => assert_eq!(what, STAPLE_RECOVERY_DOMAIN),
                other => panic!("expected the staple domain refusal, got {other:?}"),
            }
        }
        // Reversed flow has no recovery: the same opening is admitted there and is exactly the original junction.
        let reverse = junction(1.5 * limit, -500.0).unwrap();
        let original = law().characteristic_pressure(width, 1.5 * limit / width, 1.2, 1e7, 0.0, -500.0, 0.0).unwrap();
        assert_eq!(reverse.outgoing_pa, original);
        assert_eq!(reverse.recovered_pa, 0.0);
        // A coefficient or an area outside the printed range never builds.
        for (c, s) in [(-0.1, S_C), (1.1, S_C), (f64::NAN, S_C), (C_P, 0.0), (C_P, -1e-6), (C_P, f64::INFINITY)] {
            assert!(StapleRecovery::new(c, s).is_err(), "{c} {s}");
        }
    }

    #[test]
    fn g1_staple_junction_is_passive_and_its_power_balance_closes_with_the_recovery_in_it() {
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        let (width, rho) = (0.0056, 1.2);
        let (mut forward, mut reverse) = (0, 0);
        // Openings up to 0.6 mm (3.4 mm^2) stay inside S_c / sqrt(C_P) = 4.47 mm^2.
        for h in [5e-5, 1.5e-4, 3.57e-4, 6e-4] {
            for z in [3e5, 6e7, 5e9] {
                for incoming in [-300.0, 0.0, 1000.0] {
                    for mouth in [-4000.0, -50.0, 50.0, 800.0, 7400.0, 30000.0] {
                        for moving in [-8e-5, 0.0, 6e-5] {
                            let j = law().characteristic_junction_with_recovery(
                                recovery, width, h, rho, z, incoming, mouth, moving).unwrap();
                            let area = width * h;
                            let (q, v) = (j.flow_m3_s, j.flow_m3_s / area);
                            let bore = j.outgoing_pa + incoming;
                            // The characteristic junction's flow balance: wave flow = jet + swept/body flow.
                            let wave = (j.outgoing_pa - incoming) / z;
                            let flow_scale = q.abs() + moving.abs() + wave.abs();
                            assert!((q + moving - wave).abs() <= 1e-9 * flow_scale + 1e-30,
                                "h={h} z={z} in={incoming} m={mouth} mv={moving}");
                            // Independent power pieces: viscous, separated jet, and the recovery giving a part back.
                            let linear = 12.0 * 1.81e-5 * 1e-3 / h / h;
                            let viscous = linear * v * q;
                            let separation = 0.5 * rho * v * v.abs() * q;
                            let recovery_power = if q > 0.0 {
                                C_P * 0.5 * rho * (q / S_C) * (q / S_C) * q
                            } else {
                                0.0
                            };
                            let loss = viscous + separation - recovery_power;
                            let delivered = (mouth - bore) * q;
                            let scale = (mouth * q).abs() + (bore * q).abs() + viscous.abs()
                                + separation.abs() + recovery_power.abs();
                            // The power taken from the mouth and not delivered to the bore is the net loss.
                            assert!((delivered - loss).abs() <= 1e-9 * scale + 1e-30,
                                "h={h} z={z} in={incoming} m={mouth} mv={moving}: {delivered:e} vs {loss:e}");
                            // Passive: the viscous part is nonnegative, and the jet loss after the recovery is
                            // (1/2) rho v^2 q (1 - C_P S^2/S_c^2) >= 0 for forward flow.
                            assert!(loss >= -1e-9 * scale, "creates energy: {loss:e}");
                            if q > 0.0 {
                                forward += 1;
                                let net_jet = 0.5 * rho * v * v * q * (1.0 - C_P * (area / S_C) * (area / S_C));
                                assert!(net_jet > 0.0);
                                assert!((separation - recovery_power - net_jet).abs() <= 1e-9 * separation.abs());
                                assert!(j.recovered_pa > 0.0 && j.recovered_pa < j.constriction_drop_pa);
                            } else if q < 0.0 {
                                reverse += 1;
                                assert_eq!(j.recovered_pa, 0.0);
                            }
                        }
                    }
                }
            }
        }
        assert!(forward > 100 && reverse > 100, "forward {forward}, reverse {reverse}");
    }

    #[test]
    fn g0_a_shut_slit_hands_the_blade_the_whole_pressure_difference() {
        // With no flow there is no recovery: the constriction drop is the whole mouth-to-bore difference, which is what
        // presses a shut reed on its lay (found by the step-level energy test: an earlier draft returned zero here).
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        for (mouth, incoming, z, moving) in [(7400.0, 75.0, 5.9e7, 1e-7), (-300.0, 20.0, 3e5, -2e-8), (1000.0, 0.0, 6e7, 0.0)] {
            let j = law().characteristic_junction_with_recovery(
                recovery, 0.0056, 0.0, 1.2, z, incoming, mouth, moving).unwrap();
            let original = law().characteristic_pressure(0.0056, 0.0, 1.2, z, incoming, mouth, moving).unwrap();
            assert_eq!(j.outgoing_pa, original);
            assert_eq!((j.flow_m3_s, j.recovered_pa), (0.0, 0.0));
            let difference = mouth - (j.outgoing_pa + incoming);
            assert!((j.constriction_drop_pa - difference).abs() <= 1e-12 * (mouth.abs() + j.outgoing_pa.abs() + incoming.abs()),
                "{} vs {difference}", j.constriction_drop_pa);
        }
        // No difference at all: nothing presses the blade.
        // (the drive is mouth - 2 incoming - z moving = 100 - 2 * 50 - 0)
        let still = law().characteristic_junction_with_recovery(
            recovery, 0.0056, 3e-4, 1.2, 6e7, 50.0, 100.0, 0.0).unwrap();
        assert_eq!((still.flow_m3_s, still.constriction_drop_pa, still.recovered_pa), (0.0, 0.0, 0.0));
    }

    #[test]
    fn g1_vanishing_recovery_is_the_original_junction() {
        let recovery = StapleRecovery::new(1e-14, S_C).unwrap();
        let (width, rho) = (0.0056, 1.2);
        for (h, mouth, incoming, z, moving) in [
            (3.57e-4, 7400.0, 75.0, 6e7, 1e-7), (1e-4, 800.0, -30.0, 3e5, -5e-8),
            (6e-4, -2000.0, 10.0, 5e9, 0.0), (2e-4, 30000.0, 0.0, 6e7, 6e-5),
        ] {
            let j = law().characteristic_junction_with_recovery(
                recovery, width, h, rho, z, incoming, mouth, moving).unwrap();
            let original = law().characteristic_pressure(width, h, rho, z, incoming, mouth, moving).unwrap();
            assert!((j.outgoing_pa - original).abs() <= 1e-10 * (1.0 + original.abs()),
                "{} vs {original}", j.outgoing_pa);
        }
    }

    #[test]
    fn g1_condensed_viscous_recovery_retains_wall_loss_and_manufactured_flow() {
        let (width, rho, z) = (0.0056, 1.2, 6e7);
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        for opening in [1e-5, 2e-4, 6e-4] {
            for denominator in [1.0, 1.7, 100.0] {
                for velocity in [-30.0_f64, -1e-4, 0.0, 1e-4, 30.0] {
                    let flow = width * opening * velocity;
                    let wall = 12.0 * law().dynamic_viscosity_pa_s() * law().length_m()
                        * velocity / opening / opening;
                    let jet = 0.5 * rho * velocity * velocity.abs();
                    let recovered = recovery.recovered_pressure_pa(rho, flow);
                    let blade = wall + jet;
                    let total = blade - recovered;
                    // The scalar pressure force creates swept flow at the
                    // blade drop, while only the bore sees recovered head.
                    let drive = blade - recovered / denominator + (z / denominator) * flow;
                    let (actual_blade, actual_total) = law().fixed_opening_pressure_drop_with_recovery(
                        recovery, width, opening, rho, z / denominator, denominator, drive).unwrap();
                    close(actual_blade, blade, 2e-13);
                    close(actual_total, total, 2e-13);
                    let actual_flow = law().evaluate(width, opening, rho, actual_blade).unwrap();
                    close(actual_flow.flow_m3_s, flow, 2e-13);
                    assert!(actual_total * flow >= 0.0);
                    // Recovery never subtracts the no-slip linear loss.
                    assert!(actual_total * flow >= wall * flow * (1.0 - 2e-13));
                }
            }
        }
    }

    #[test]
    fn g0_condensed_viscous_recovery_keeps_identity_limits_and_refusals() {
        let recovery = StapleRecovery::new(C_P, S_C).unwrap();
        let zero = StapleRecovery::new(0.0, S_C).unwrap();
        for (opening, drive) in [(2e-4, -100.0), (0.0, 2000.0), (2e-4, 0.0)] {
            let old = law().fixed_opening_pressure_drop(0.0056, opening, 1.2, 1e7, drive).unwrap();
            for selected in [recovery, zero] {
                let actual = law().fixed_opening_pressure_drop_with_recovery(
                    selected, 0.0056, opening, 1.2, 1e7, 2.0, drive).unwrap();
                assert_eq!((actual.0.to_bits(), actual.1.to_bits()), (old.to_bits(), old.to_bits()));
            }
        }
        let old = law().fixed_opening_pressure_drop(0.0056, 2e-4, 1.2, 1e7, 500.0).unwrap();
        let actual = law().fixed_opening_pressure_drop_with_recovery(
            zero, 0.0056, 2e-4, 1.2, 1e7, 2.0, 500.0).unwrap();
        assert_eq!((actual.0.to_bits(), actual.1.to_bits()), (old.to_bits(), old.to_bits()));
        // Mechanical condensation cannot grant a nonpassive diffuser domain.
        assert!(matches!(law().fixed_opening_pressure_drop_with_recovery(
            recovery, 0.02, 1e-3, 1.2, 1e7, 100.0, 500.0),
            Err(AcousticRealizeError::Reed { what: STAPLE_RECOVERY_DOMAIN })));
        for denominator in [0.0, 0.5, f64::NAN, f64::INFINITY] {
            assert!(law().fixed_opening_pressure_drop_with_recovery(
                recovery, 0.0056, 2e-4, 1.2, 1e7, denominator, 500.0).is_err());
        }
    }
}
