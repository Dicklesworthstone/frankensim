//! A closed, evolving crystal inventory against one constant pressure load.
//!
//! The original caloric/density law determines actual parcel volume. Reference
//! volume remains a separate immutable mechanical quadrature convention. The
//! parcel advances the shared phase history and accepts its energy/work closure
//! atomically; it supplies no moving mesh, piston inertia or pressure EOS.

use super::super::{
    CrystallizationAdaptivity, CrystallizationError, CrystallizationLaw, CrystallizationState,
    EquivalentAgeRemelting,
};
use super::{
    IsobaricCrystallization, IsobaricCrystallizationControls, IsobaricCrystallizationStep,
};

/// Accepted mass, phase history and actual volume at a fixed external pressure.
/// Clones preserve the original immutable sources and the complete retry state.
#[derive(Clone)]
pub struct IsobaricCrystallizationParcel<'law, 'source> {
    closure: IsobaricCrystallization<'law, 'source>,
    state: CrystallizationState,
    mass_kg: f64,
    reference_volume_m3: f64,
    volume_m3: f64,
    time_s: f64,
    remelting: Option<EquivalentAgeRemelting>,
}

impl<'law, 'source> IsobaricCrystallizationParcel<'law, 'source> {
    /// Adopt an accepted kinetic state, conserved mass and explicit mechanical
    /// reference volume. The local elapsed-time clock starts at this adoption;
    /// the supplied exposure is retained, including any earlier phase history.
    ///
    /// The actual initial volume is m*v(T,X), independently of the supplied
    /// reference volume. The original caloric, density and kinetic sources must
    /// all apply at this constant absolute pressure and composition. This is a
    /// quasistatic pressure-load model, not a rigid sealed-volume caloric law.
    pub fn new(
        law: &'law CrystallizationLaw<'source>,
        state: CrystallizationState,
        mass_kg: f64,
        reference_volume_m3: f64,
        pressure_pa: f64,
        controls: IsobaricCrystallizationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, CrystallizationError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if !mass_kg.is_finite()
            || mass_kg <= 0.0
            || !reference_volume_m3.is_finite()
            || reference_volume_m3 <= 0.0
            || !(reference_volume_m3 / mass_kg).is_finite()
            || reference_volume_m3 / mass_kg <= 0.0
        {
            return Err(CrystallizationError::Invalid(
                "isobaric parcel requires representable positive mass and reference volume",
            ));
        }
        let closure = IsobaricCrystallization::new(law, pressure_pa)?;
        let volume_m3 = closure.volume_m3(state, mass_kg, controls.caloric, &mut cancelled)?;
        poll(&mut cancelled)?;
        Ok(Self {
            closure,
            state,
            mass_kg,
            reference_volume_m3,
            volume_m3,
            time_s: 0.0,
            remelting: None,
        })
    }

    /// Select the original equivalent-age remelting model once. Growth alone
    /// remains the default; selection does not alter the accepted physical state.
    pub fn with_equilibrium_remelting(
        mut self,
        model: EquivalentAgeRemelting,
    ) -> Result<Self, CrystallizationError> {
        if self.remelting.is_some() {
            return Err(CrystallizationError::Invalid(
                "isobaric parcel remelting requires a unique selection",
            ));
        }
        if self.law().equivalent_age_remelting()? != model {
            return Err(CrystallizationError::StateMismatch);
        }
        self.remelting = Some(model);
        Ok(self)
    }

    pub fn law(&self) -> &'law CrystallizationLaw<'source> {
        self.closure.law()
    }
    pub const fn state(&self) -> CrystallizationState {
        self.state
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    /// Fixed mechanical reference volume; never replaced by free phase volume.
    pub const fn reference_volume_m3(&self) -> f64 {
        self.reference_volume_m3
    }
    pub fn pressure_pa(&self) -> f64 {
        self.closure.pressure_pa()
    }
    /// Actual current volume from the original densities and kinetic fraction.
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    /// Physical elapsed time since adoption of the initial accepted state.
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn remelting_model(&self) -> Option<EquivalentAgeRemelting> {
        self.remelting
    }

    /// Total enthalpy on the source datum. Increment reports remain the accurate
    /// account when a large arbitrary datum hides a small physical heat change.
    pub fn enthalpy_j(&self) -> Result<f64, CrystallizationError> {
        finite(self.mass_kg * self.state.specific_enthalpy_j_kg())
    }

    /// Total U=H-pV on the declared pressure/reference convention. This is not
    /// an independent variable-pressure EOS. Use step increments for balances.
    pub fn internal_energy_j(&self) -> Result<f64, CrystallizationError> {
        finite(self.enthalpy_j()? - finite(self.pressure_pa() * self.volume_m3)?)
    }

    /// Propose signed heat and its quasistatic volume/work response at unchanged
    /// kinetic time. No finite-speed boundary trajectory or momentum is inferred.
    /// Positive work leaves the parcel for its declared pressure-load owner.
    pub fn apply_heat(
        &self,
        heat_j: f64,
        controls: IsobaricCrystallizationControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, IsobaricCrystallizationStep), CrystallizationError> {
        poll(&mut cancelled)?;
        let step = self.closure.apply_heat(
            self.state,
            self.mass_kg,
            heat_j,
            controls,
            self.remelting,
            &mut cancelled,
        )?;
        let next = self.candidate(&step, self.time_s)?;
        poll(&mut cancelled)?;
        Ok((next, step))
    }

    /// Advance the actual phase history over a positive elapsed interval while
    /// adding signed heat linearly in time. Adiabatic growth is supported: zero
    /// heat does not freeze kinetics, and contraction can return pressure work.
    ///
    /// The shared adaptive kernel owns all time integration and its trial caps.
    /// The report's duration matches the represented increment of this clock.
    /// Source, phase, work, accuracy or final cancellation failure leaves this
    /// entire accepted parcel unchanged, including volume and elapsed time.
    pub fn advance_heat(
        &self,
        heat_j: f64,
        duration_s: f64,
        controls: IsobaricCrystallizationControls,
        adaptivity: CrystallizationAdaptivity,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(Self, IsobaricCrystallizationStep), CrystallizationError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(CrystallizationError::Invalid(
                "isobaric parcel advance requires positive elapsed time",
            ));
        }
        let time_s = finite(self.time_s + duration_s)?;
        let represented_duration = time_s - self.time_s;
        if represented_duration <= 0.0 {
            return Err(CrystallizationError::Invalid(
                "isobaric parcel elapsed-time increment is unrepresentable",
            ));
        }
        let step = self.closure.advance_heat(
            self.state,
            self.mass_kg,
            heat_j,
            represented_duration,
            controls,
            adaptivity,
            self.remelting,
            &mut cancelled,
        )?;
        let next = self.candidate(&step, time_s)?;
        poll(&mut cancelled)?;
        Ok((next, step))
    }

    fn candidate(
        &self,
        step: &IsobaricCrystallizationStep,
        time_s: f64,
    ) -> Result<Self, CrystallizationError> {
        // This exact check rejects a mutable density source after admission;
        // a tolerance must not turn changed constitutive data into state drift.
        if step.initial_volume_m3 != self.volume_m3 {
            return Err(CrystallizationError::Invalid(
                "original density source changed after isobaric parcel acceptance",
            ));
        }
        let mut next = self.clone();
        next.state = step.state;
        next.volume_m3 = step.final_volume_m3;
        next.time_s = time_s;
        Ok(next)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CrystallizationError> {
    if cancelled() {
        Err(CrystallizationError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: f64) -> Result<f64, CrystallizationError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CrystallizationError::Invalid(
            "isobaric parcel total is unrepresentable",
        ))
    }
}

#[cfg(test)]
#[path = "parcel/tests.rs"]
mod tests;
