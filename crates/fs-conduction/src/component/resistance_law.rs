//! A contact resistance that follows the temperatures of its two traces.
//!
//! A [`TraceContact`](super::TraceContact) states the owner's contact between
//! two trace ports at one resistance. [`TraceResistanceLaw`] is that
//! resistance as the card gives it at the temperatures of the two traces: at
//! every quadrature point the contact's conductance is its conductance at the
//! law's reference resistance times [`TraceResistanceLaw::factor`], which is
//! `R''_ref / R''(T_a, T_b)`. fs-system's `ConductanceMatrix::following`
//! takes that factor through its `ConductanceLaw` hook, and the element then
//! is the owner's pointwise temperature contact
//! (`interface::temperature::TemperatureContact::new_pointwise` and
//! `new_nonmatching`) between two bodies on the kernel.
//!
//! The card is asked as the owner asks it. A claim is selected once, at the
//! reference point, and pinned: a later state cannot switch claims where the
//! claim's support ends. Every other coordinate of the reference point
//! (pressure, finish, medium) stays where it was declared. A state the pinned
//! claim has no answer for is refused with fs-matdb's reason; nothing is
//! clamped or extrapolated.
//!
//! The slope is the one the owner's Newton tangent takes
//! (`TemperatureContacts::linearize_temperature_for_newton`): of a linear
//! curve, the right segment's at an interior knot and the inward segment's
//! at either end, from the same function.
//!
//! No claim: a claim that is a surface over a temperature coordinate, and
//! any interpolation other than a constant scalar or a linear curve, is
//! refused when the law is bound. No material state (a liquid fraction)
//! enters the query.

use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyValue, QueryPoint, SelectionPolicy,
};

use crate::ConductionError;
use crate::interface::InterfaceResistance;
use crate::interface::temperature::TemperatureContactCoordinates;
use crate::interface::temperature::tangent::newton_segment_slope;

/// How the pinned claim's resistance moves with the two trace temperatures.
#[derive(Debug, Clone, PartialEq)]
enum Slope {
    /// It does not: a scalar, or a claim over coordinates that are not the
    /// law's temperatures.
    Level,
    /// A linear curve over the temperature coordinate `axis`, and the share
    /// of either trace's temperature in that coordinate.
    Curve {
        axis: String,
        knots: Vec<(f64, f64)>,
        shares: [f64; 2],
    },
}

/// An interface card's area-specific contact resistance at the temperatures
/// of two traces, relative to its value at a reference point.
#[derive(Debug, Clone, PartialEq)]
pub struct TraceResistanceLaw {
    name: String,
    card: InterfaceSystemCard,
    reference_point: QueryPoint,
    coordinates: TemperatureContactCoordinates,
    reference: InterfaceResistance,
    slope: Slope,
}

impl TraceResistanceLaw {
    /// Selects one claim of `card` at `point` under `selection` and pins it.
    ///
    /// `coordinates` names the card's temperature coordinates: one axis that
    /// follows the mean of the two trace temperatures, or two that follow
    /// side A's and side B's. `point` holds a reference value of each and
    /// every other coordinate the card asks for. The contact itself is built
    /// at [`Self::reference_resistance`].
    ///
    /// # Errors
    /// The owner's refusals of a temperature contact's reference: a blank or
    /// repeated temperature axis, a reference point without a positive
    /// absolute temperature on each, and fs-matdb's refusal of the query. A
    /// claim this law does not carry (module documentation) is refused too.
    pub fn bind(
        name: impl Into<String>,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
        coordinates: TemperatureContactCoordinates,
        selection: SelectionPolicy,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        coordinates.admit_reference(&name, point)?;
        let reference = InterfaceResistance::from_card(&name, card, point, selection)?;
        let claim = card
            .claims()
            .claim(reference.receipt().selected)
            .ok_or_else(|| refused(&name, "the selected contact claim is missing"))?;
        let shares = |axis: &str| -> [f64; 2] {
            match &coordinates {
                TemperatureContactCoordinates::Mean(temperature) => {
                    if axis == temperature {
                        [0.5; 2]
                    } else {
                        [0.0; 2]
                    }
                }
                TemperatureContactCoordinates::TwoSided(axes) => axes
                    .each_ref()
                    .map(|temperature| if axis == temperature { 1.0 } else { 0.0 }),
            }
        };
        // Whether a claim's axis is one of the law's temperature coordinates.
        let follows = |axis: &String| coordinates.axes().contains(axis);
        let slope = match (&claim.value, &claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {
                Slope::Level
            }
            (
                PropertyValue::Curve {
                    abscissa, knots, ..
                },
                InterpolationPolicy::LinearInside,
            ) => {
                if follows(abscissa) {
                    Slope::Curve {
                        axis: abscissa.clone(),
                        knots: knots.clone(),
                        shares: shares(abscissa),
                    }
                } else {
                    Slope::Level
                }
            }
            (PropertyValue::Surface2D { axes, .. }, InterpolationPolicy::LinearInside)
                if !axes.iter().any(follows) =>
            {
                Slope::Level
            }
            _ => {
                return Err(refused(
                    &name,
                    "a trace resistance law carries a constant scalar or a linear curve; a surface over a temperature coordinate, and any other interpolation, is not built",
                ));
            }
        };
        Ok(Self {
            name,
            card: card.clone(),
            reference_point: point.clone(),
            coordinates,
            reference,
            slope,
        })
    }

    /// The interface's name, as the card's refusals print it.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The card's answer at the reference point, with its receipt: the
    /// resistance the contact is built at, where [`Self::factor`] is one.
    #[must_use]
    pub const fn reference_resistance(&self) -> &InterfaceResistance {
        &self.reference
    }

    /// The card's temperature coordinates.
    #[must_use]
    pub const fn coordinates(&self) -> &TemperatureContactCoordinates {
        &self.coordinates
    }

    /// `R''_ref / R''(T_a, T_b)`, the factor on a contact point's conductance
    /// at the temperatures `kelvin_a` and `kelvin_b` of the two traces there,
    /// and its partial derivatives with respect to them [1/K].
    ///
    /// # Errors
    /// [`ConductionError::Interface`] for a temperature that is not positive
    /// and finite, and for a state the pinned claim has no answer for, with
    /// the reason fs-matdb gave.
    pub fn factor(&self, kelvin_a: f64, kelvin_b: f64) -> Result<(f64, [f64; 2]), ConductionError> {
        let point = self.coordinates.query_point_at(
            &self.name,
            &self.reference_point,
            [kelvin_a, kelvin_b],
        )?;
        let here = InterfaceResistance::from_card_pinned(
            &self.name,
            &self.card,
            &point,
            self.reference.receipt().selected,
        )?;
        let resistance = here.value_m2_k_per_w();
        let factor = self.reference.value_m2_k_per_w() / resistance;
        let partials = match &self.slope {
            Slope::Level => [0.0; 2],
            Slope::Curve {
                axis,
                knots,
                shares,
            } => {
                let value = *point.axes().get(axis).ok_or_else(|| {
                    refused(
                        &self.name,
                        "the law's curve lost its temperature coordinate",
                    )
                })?;
                let slope = newton_segment_slope(&self.name, knots, value)?;
                // d(R_ref / R) = -(R_ref / R) dR / R.
                shares.map(|share| -(share * slope / resistance) * factor)
            }
        };
        Ok((factor, partials))
    }

    /// Words that identify the law: a hash of the card, the pinned claim, the
    /// temperature coordinates and the reference point, as four words.
    pub fn identity_words(&self, words: &mut Vec<u64>) {
        fn text(hasher: &mut fs_blake3::DomainHasher, text: &str) {
            hasher.update(&(text.len() as u64).to_le_bytes());
            hasher.update(text.as_bytes());
        }
        let mut hasher = fs_blake3::DomainHasher::new("fs-conduction/trace-resistance-law/v1");
        hasher.update(self.card.content_hash().as_bytes());
        hasher.update(self.reference.receipt().selected.0.as_bytes());
        let axes = self.coordinates.axes();
        hasher.update(&(axes.len() as u64).to_le_bytes());
        for axis in axes {
            text(&mut hasher, axis);
        }
        for (axis, value) in self.reference_point.axes() {
            text(&mut hasher, axis);
            hasher.update(&value.to_bits().to_le_bytes());
        }
        let hash = hasher.finalize();
        for word in hash.as_bytes().as_chunks::<8>().0 {
            words.push(u64::from_le_bytes(*word));
        }
    }
}

fn refused(name: &str, what: &str) -> ConductionError {
    ConductionError::Interface {
        interface: name.into(),
        what: what.into(),
        fix: "bind the law to a claim that is a constant scalar or a linear curve over the declared temperature coordinate, and keep both traces inside it".into(),
    }
}
