//! Conservative wall-opening exchange between two hydrostatic, two-layer columns.
//!
//! The opening is split at both layer interfaces and every neutral-pressure plane.
//! Each slab integrates the existing quadratic orifice law analytically using its
//! donor density. Opposite transfers remain separate even when their net mass
//! flow vanishes. No compartment state, species, or energy is advanced here.
//!
//! Reference: Peacock et al., NIST TN 1889v1 (2015), section 4.1, equations
//! 4.1–4.4, <https://doi.org/10.6028/NIST.TN.1889v1>. The caller supplies Cd;
//! the reference's coefficient is not substituted for the original airflow law.
//! This is a low-Mach, piecewise-constant-density engineering closure. It does
//! not model compressible choking, horizontal ceiling vents, wind fields, jet
//! entrainment, or empirical validation. Retained coefficient uncertainties and
//! validity cards are declarations, not propagated uncertainty or admission.

use crate::{
    AirflowError, LossElement, LossResistance, SHARP_EDGED_ORIFICE_CD, SourceProvenance,
    ToleranceBasis, sharp_edged_orifice_loss,
};
use fs_math::det;
use fs_qty::{
    Acceleration, Area, Density, Dimensionless, Length, MassFlowRate, Power, Pressure, Qty,
};

mod flow;
pub use flow::wall_vent_flow;

/// At most three constant-density intervals, each with at most one neutral plane.
pub const MAX_VENT_SLABS: usize = 6;

/// Specific enthalpy in J/kg; the caller declares the common energy reference.
pub type SpecificEnthalpy = Qty<2, 0, -2, 0, 0, 0>;

/// Layer identity in a column, independent of the direction of transfer.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Layer {
    /// Gas below the interface.
    #[default]
    Lower,
    /// Gas above the interface.
    Upper,
}
impl Layer {
    /// Canonical index into `[lower, upper]` arrays.
    pub const fn index(self) -> usize {
        match self {
            Self::Lower => 0,
            Self::Upper => 1,
        }
    }
}

/// One admitted hydrostatic column. Heights use a shared world vertical datum.
///
/// Pressure is absolute at `reference_height`; density is constant within each
/// layer. The interface may lie outside the opening, including an absent layer.
/// Equal densities represent an unstratified atmosphere without a special case.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoLayerColumn {
    reference_height: Length,
    reference_pressure: Pressure,
    interface_height: Length,
    densities: [Density; 2],
}
impl TwoLayerColumn {
    /// Bind finite heights, positive absolute pressure, and positive densities.
    pub fn new(
        reference_height: Length,
        reference_pressure: Pressure,
        interface_height: Length,
        lower_density: Density,
        upper_density: Density,
    ) -> Result<Self, VentError> {
        if !reference_height.value().is_finite() || !interface_height.value().is_finite() {
            return Err(VentError::Invalid(
                "finite reference and interface heights required",
            ));
        }
        if !positive(reference_pressure.value())
            || !positive(lower_density.value())
            || !positive(upper_density.value())
        {
            return Err(VentError::Invalid(
                "positive absolute pressure and layer densities required",
            ));
        }
        Ok(Self {
            reference_height,
            reference_pressure,
            interface_height,
            densities: [lower_density, upper_density],
        })
    }

    /// Height at which the absolute reference pressure is prescribed.
    pub const fn reference_height(self) -> Length {
        self.reference_height
    }
    /// Retained absolute reference pressure.
    pub const fn reference_pressure(self) -> Pressure {
        self.reference_pressure
    }
    /// Height separating the lower and upper layers.
    pub const fn interface_height(self) -> Length {
        self.interface_height
    }
    /// Retained density of one layer.
    pub const fn density(self, layer: Layer) -> Density {
        self.densities[layer.index()]
    }

    /// Evaluate absolute hydrostatic pressure, refusing a nonpositive result.
    /// Gravity is the nonnegative magnitude of downward acceleration.
    pub fn pressure_at(self, height: Length, gravity: Acceleration) -> Result<Pressure, VentError> {
        validate_gravity(gravity)?;
        if !height.value().is_finite() {
            return Err(VentError::Invalid("finite evaluation height required"));
        }
        if gravity.value() == 0.0 {
            return Ok(self.reference_pressure);
        }
        let pressure =
            self.reference_pressure.value() - gravity.value() * self.density_height(height)?;
        if !positive(pressure) {
            return Err(VentError::Invalid(
                "hydrostatic column left positive absolute pressure",
            ));
        }
        Ok(Pressure::new(pressure))
    }

    fn density_height(self, height: Length) -> Result<f64, VentError> {
        let z = height.value();
        if !z.is_finite() {
            return Err(VentError::Invalid("finite evaluation height required"));
        }
        let reference = self.reference_height.value();
        let interface = self.interface_height.value();
        let lower = z.min(interface) - reference.min(interface);
        let upper = z.max(interface) - reference.max(interface);
        finite(
            self.densities[0].value() * lower + self.densities[1].value() * upper,
            "hydrostatic density-height integral",
        )
    }
}

/// An explicitly sourced discharge coefficient lowered into the original loss law.
///
/// The retained loss is referenced to area 1 m² and density 1 kg/m³, so
/// `R0 = 1/(2 Cd²)`. Actual geometry and donor density scale this same resistance;
/// no differently normalized Cd or independent flow law is introduced.
#[derive(Clone, Debug, PartialEq)]
pub struct DischargeCoefficient {
    value: Dimensionless,
    reference_loss: LossElement,
}
impl DischargeCoefficient {
    /// Admit `0 < Cd <= 1` and explicit authority for resistance uncertainty.
    /// The uncertainty applies to resistance, as in [`LossElement`], not to Cd.
    pub fn new(
        value: Dimensionless,
        source: SourceProvenance,
        relative_resistance_uncertainty: f64,
        basis: ToleranceBasis,
    ) -> Result<Self, VentError> {
        if !positive(value.value()) || value.value() > 1.0 {
            return Err(VentError::Invalid(
                "discharge coefficient must lie in (0,1]",
            ));
        }
        if source.citation.trim().is_empty() || source.identifier.trim().is_empty() {
            return Err(VentError::Invalid(
                "discharge coefficient requires explicit source identity",
            ));
        }
        let resistance = 0.5 / (value.value() * value.value());
        let reference_loss = LossElement::new(
            "stratified opening",
            LossResistance::new(resistance),
            relative_resistance_uncertainty,
            source,
            basis,
        )?;
        Ok(Self {
            value,
            reference_loss,
        })
    }

    /// Retain the crate's original sharp-edged coefficient, source, uncertainty,
    /// and Reynolds validity declaration by calling its existing lowering path.
    pub fn sharp_edged() -> Result<Self, VentError> {
        Ok(Self {
            value: Dimensionless::new(SHARP_EDGED_ORIFICE_CD),
            reference_loss: sharp_edged_orifice_loss(
                "stratified sharp-edged opening",
                Area::new(1.0),
                Density::new(1.0),
            )?,
        })
    }
    /// Declared coefficient, with its original dimensionless normalization.
    pub const fn value(&self) -> Dimensionless {
        self.value
    }
    /// Original quadratic loss and retained model authority at unit area/density.
    pub const fn reference_loss(&self) -> &LossElement {
        &self.reference_loss
    }
}

/// A rectangular vertical opening; zero width or height explicitly seals it.
#[derive(Clone, Debug, PartialEq)]
pub struct WallOpening {
    bottom: Length,
    top: Length,
    width: Length,
    discharge: DischargeCoefficient,
}
impl WallOpening {
    /// Bind finite geometry and the retained discharge law. Negative width or
    /// reversed vertical bounds refuse. Heights use both columns' world datum.
    pub fn new(
        bottom: Length,
        top: Length,
        width: Length,
        discharge: DischargeCoefficient,
    ) -> Result<Self, VentError> {
        let height = top.value() - bottom.value();
        if !bottom.value().is_finite()
            || !top.value().is_finite()
            || !height.is_finite()
            || height < 0.0
            || !width.value().is_finite()
            || width.value() < 0.0
        {
            return Err(VentError::Invalid(
                "finite ordered opening heights and nonnegative width required",
            ));
        }
        let area = width.value() * height;
        if !area.is_finite() || (area == 0.0 && width.value() != 0.0 && height != 0.0) {
            return Err(VentError::Invalid("opening area is unrepresentable"));
        }
        Ok(Self {
            bottom,
            top,
            width,
            discharge,
        })
    }
    /// Bottom edge above the common vertical datum.
    pub const fn bottom(&self) -> Length {
        self.bottom
    }
    /// Top edge above the common vertical datum.
    pub const fn top(&self) -> Length {
        self.top
    }
    /// Horizontal open width, including any explicitly prescribed closure.
    pub const fn width(&self) -> Length {
        self.width
    }
    /// Geometric open area.
    pub fn area(&self) -> Area {
        Area::new(self.width.value() * (self.top.value() - self.bottom.value()))
    }
    /// Original source-bearing discharge declaration.
    pub const fn discharge(&self) -> &DischargeCoefficient {
        &self.discharge
    }
}

/// One nonzero, single-direction slab between interfaces or neutral planes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabFlow {
    /// Lower integration height [m].
    pub bottom: Length,
    /// Upper integration height [m].
    pub top: Length,
    /// Integrated vertical extent [m], computed in local coordinates. Near a
    /// neutral plane, subtracting rounded world-coordinate endpoints can lose
    /// this small positive extent; the flow integral uses this value directly.
    pub height: Length,
    /// Layer intersected on the left, regardless of direction.
    pub left_layer: Layer,
    /// Layer intersected on the right, regardless of direction.
    pub right_layer: Layer,
    /// Signed mass rate [kg/s], positive left to right.
    pub mass_rate: MassFlowRate,
    /// Density of the actual donor layer [kg/m³].
    pub upstream_density: Density,
    /// Left-minus-right pressure at the lower endpoint [Pa].
    pub pressure_difference_bottom: Pressure,
    /// Left-minus-right pressure at the upper endpoint [Pa].
    pub pressure_difference_top: Pressure,
}
impl SlabFlow {
    const EMPTY: Self = Self {
        bottom: Length::new(0.0),
        top: Length::new(0.0),
        height: Length::new(0.0),
        left_layer: Layer::Lower,
        right_layer: Layer::Lower,
        mass_rate: MassFlowRate::new(0.0),
        upstream_density: Density::new(0.0),
        pressure_difference_bottom: Pressure::new(0.0),
        pressure_difference_top: Pressure::new(0.0),
    };
}

/// Per-layer mass source rates; positive means entering that layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerMassRates {
    /// Left column, ordered `[lower, upper]` [kg/s].
    pub left: [MassFlowRate; 2],
    /// Right column, ordered `[lower, upper]` [kg/s].
    pub right: [MassFlowRate; 2],
}

/// Per-layer advected enthalpy rates; positive means entering that layer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LayerEnthalpyRates {
    /// Left column, ordered `[lower, upper]` [W].
    pub left: [Power; 2],
    /// Right column, ordered `[lower, upper]` [W].
    pub right: [Power; 2],
}

/// Bounded, allocation-free exchange result, in increasing slab-height order.
#[derive(Clone, Debug, PartialEq)]
pub struct VentFlow {
    slabs: [SlabFlow; MAX_VENT_SLABS],
    count: usize,
    rates: LayerMassRates,
}
impl VentFlow {
    /// Every nonzero transfer, retaining both directions and both layer labels.
    pub fn slabs(&self) -> &[SlabFlow] {
        &self.slabs[..self.count]
    }
    /// Equal-and-opposite transfer accounting aggregated by receiving layer.
    pub const fn mass_rates(&self) -> LayerMassRates {
        self.rates
    }

    /// Advect one common-reference donor enthalpy per layer without evaluating
    /// or replacing the caller's thermodynamic model. Sensible or total enthalpy
    /// is accepted; all four values must be finite and use the same energy gauge.
    /// Each slab's donor power is added to one layer and subtracted from the
    /// other. Net mass flow must not be used as a substitute in counterflow.
    pub fn advect_enthalpy(
        &self,
        left: [SpecificEnthalpy; 2],
        right: [SpecificEnthalpy; 2],
    ) -> Result<LayerEnthalpyRates, VentError> {
        if !left.iter().chain(&right).all(|h| h.value().is_finite()) {
            return Err(VentError::Invalid(
                "finite donor specific enthalpies required",
            ));
        }
        let mut rates = [Sum::default(); 4];
        for slab in self.slabs() {
            let l = slab.left_layer.index();
            let r = slab.right_layer.index();
            let mass = slab.mass_rate.value();
            let h = if mass > 0.0 { left[l] } else { right[r] };
            let power = finite(mass * h.value(), "donor enthalpy rate")?;
            rates[l].add(-power)?;
            rates[2 + r].add(power)?;
        }
        Ok(LayerEnthalpyRates {
            left: [Power::new(rates[0].value), Power::new(rates[1].value)],
            right: [Power::new(rates[2].value), Power::new(rates[3].value)],
        })
    }
}

/// Invalid data, unrepresentable arithmetic, original loss refusal, or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum VentError {
    /// Input or physical domain is not admitted.
    Invalid(&'static str),
    /// An actual intermediate or aggregate cannot be represented finitely.
    NonFinite(&'static str),
    /// The original quadratic-loss owner refused the declaration.
    Loss(AirflowError),
    /// Cancellation occurred before a complete result was published.
    Cancelled,
}
impl core::fmt::Display for VentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "stratified wall opening: {self:?}")
    }
}
impl std::error::Error for VentError {}
impl From<AirflowError> for VentError {
    fn from(error: AirflowError) -> Self {
        Self::Loss(error)
    }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn finite(value: f64, stage: &'static str) -> Result<f64, VentError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(VentError::NonFinite(stage))
    }
}
fn validate_gravity(gravity: Acceleration) -> Result<(), VentError> {
    if !gravity.value().is_finite() || gravity.value() < 0.0 {
        Err(VentError::Invalid(
            "finite nonnegative downward gravity required",
        ))
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), VentError> {
    if cancelled() {
        Err(VentError::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Clone, Copy, Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), VentError> {
        let y = finite(term - self.correction, "layer transfer summand")?;
        let next = finite(self.value + y, "layer transfer sum")?;
        self.correction = finite((next - self.value) - y, "layer transfer compensation")?;
        self.value = next;
        Ok(())
    }
}

#[cfg(test)]
mod tests;
