//! A second finite source parcel on the opposing side of the same membrane.

use super::*;

/// Explicit opposing chamber at the same membrane reference shape.
///
/// Its volume decreases by exactly the sweep which increases the first parcel's
/// volume. The remaining walls are rigid, and no mass or heat crosses the sheet.
/// Independent entropy controls admit different source-temperature domains.
#[derive(Clone, Debug, PartialEq)]
pub struct SealedMembraneExterior {
    pub gas: MembraneGas,
    pub volume_interval_m3: [f64; 2],
    pub entropy: GasEntropySolve,
}

/// Same-shape opposing source state and its reference-relative native account.
#[derive(Clone, Debug, PartialEq)]
pub struct ExteriorGasEvaluation {
    pub gas: MembraneGas,
    pub gas_volume: GasVolumeStep,
    pub pressure: GasPressureResponse,
}

pub(super) enum Exterior {
    Pressure(f64),
    Sealed(SealedMembraneExterior),
}

impl Exterior {
    pub(super) fn trial(
        &self,
        swept_m3: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Option<ExteriorGasEvaluation>, GasMembraneError> {
        let Self::Sealed(chamber) = self else {
            return Ok(None);
        };
        let reference = chamber.gas.volume_m3();
        let volume = reference - swept_m3;
        if !volume.is_finite()
            || volume < chamber.volume_interval_m3[0]
            || volume > chamber.volume_interval_m3[1]
            || (swept_m3 != 0.0 && volume == reference)
        {
            return Err(GasMembraneError::Invalid(
                "opposing swept volume is outside its interval or unrepresentable",
            ));
        }
        let (gas, gas_volume) = chamber.gas.at_volume(volume, chamber.entropy, cancelled)?;
        let pressure = gas.pressure_response()?;
        Ok(Some(ExteriorGasEvaluation {
            gas,
            gas_volume,
            pressure,
        }))
    }

    pub(super) fn pressure_pa(&self, gas: Option<&ExteriorGasEvaluation>) -> f64 {
        match self {
            Self::Pressure(pressure) => *pressure,
            Self::Sealed(_) => gas.expect("same-shape exterior trial").pressure.pressure_pa,
        }
    }

    pub(super) fn reservoir_work_j(&self, swept_m3: f64) -> f64 {
        match self {
            Self::Pressure(pressure) => pressure * swept_m3,
            Self::Sealed(_) => 0.0,
        }
    }

    pub(super) fn admit_modes(
        &self,
        value: Option<&ExteriorGasEvaluation>,
    ) -> Result<(), GasMembraneError> {
        let Self::Sealed(chamber) = self else {
            return Ok(());
        };
        let gas = &value.expect("accepted exterior gas").gas;
        let volume = gas.volume_m3();
        let temperature = match gas {
            MembraneGas::Ideal(gas) => gas.state().state().temperature,
            MembraneGas::PengRobinson(gas) => gas.state().temperature_k(),
        };
        if volume <= chamber.volume_interval_m3[0]
            || volume >= chamber.volume_interval_m3[1]
            || temperature <= chamber.entropy.minimum_temperature_k
            || temperature >= chamber.entropy.maximum_temperature_k
        {
            return Err(GasMembraneError::Invalid(
                "opposing gas modes need equilibrium strictly inside its independent volume and temperature brackets",
            ));
        }
        Ok(())
    }
}

impl<'a> GasMembrane<'a> {
    /// Bind two original closed parcels to opposite sides of one membrane.
    ///
    /// Positive orientation points out of the first gas into `exterior`:
    /// V1 = V1_reference + sweep, V2 = V2_reference - sweep. Both parcel states
    /// use that SAME geometric sweep; no evolving-volume drift is accumulated.
    /// Each chamber retains its own source, mass, entropy and positive volume
    /// interval. The second parcel uses `exterior.entropy`; existing evaluation,
    /// equilibrium and motion controls continue to govern the first parcel.
    ///
    /// This replaces the infinite fixed-pressure reservoir. All existing solve,
    /// mode and motion methods use both gases through the common force/tangent.
    #[allow(clippy::too_many_arguments)]
    pub fn with_sealed_exterior(
        patch: &'a MembranePatch,
        gas: impl Into<MembraneGas>,
        reference_m: &[[f64; 3]],
        fixed: &[Option<f64>],
        volume_interval_m3: [f64; 2],
        exterior: SealedMembraneExterior,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, GasMembraneError> {
        let mut model = Self::new(
            patch,
            gas,
            reference_m,
            fixed,
            0.0,
            volume_interval_m3,
            &mut cancelled,
        )?;
        let [lo, hi] = exterior.volume_interval_m3;
        if !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || lo >= hi
            || exterior.gas.volume_m3() < lo
            || exterior.gas.volume_m3() > hi
            || !(model.gas.volume_m3() + exterior.gas.volume_m3()).is_finite()
        {
            return Err(GasMembraneError::Invalid(
                "positive opposing volume interval containing its gas and representable total volume required",
            ));
        }
        // The original source owner validates its independent entropy controls
        // and reference temperature; no duplicate bracket or caloric law.
        let boundary = Exterior::Sealed(exterior);
        boundary.trial(0.0, &mut cancelled)?;
        poll(&mut cancelled)?;
        model.exterior = boundary;
        Ok(model)
    }

    /// Immutable opposing inventory, geometric interval and entropy controls.
    /// `None` means the original held-pressure exterior boundary.
    #[must_use]
    pub fn sealed_exterior(&self) -> Option<&SealedMembraneExterior> {
        match &self.exterior {
            Exterior::Pressure(_) => None,
            Exterior::Sealed(chamber) => Some(chamber),
        }
    }
}

impl ExteriorGasEvaluation {
    pub(super) fn add_bulk_tangent(
        &self,
        gradient: &[f64],
        direction: &[f64],
        action: &mut [f64],
    ) -> Result<(), GasMembraneError> {
        let opposing_gradient: Vec<_> = gradient.iter().map(|g| -g).collect();
        let zero_hessian = vec![0.0; gradient.len()];
        // The first parcel's native action already contains (p1-p2) Hess(V1).
        // Pull back only the second parcel's ORIGINAL compressibility: its
        // gradient is -grad(V1), and its local pressure jump is zero here.
        // This avoids counting exterior geometric traction twice or subtracting
        // two large absolute-pressure Hessian actions at near-zero net load.
        let bulk = self.gas.tangent(
            &opposing_gradient,
            self.pressure.pressure_pa,
            direction,
            &zero_hessian,
        )?;
        for (value, extra) in action.iter_mut().zip(bulk) {
            *value += extra;
        }
        Ok(())
    }
}
