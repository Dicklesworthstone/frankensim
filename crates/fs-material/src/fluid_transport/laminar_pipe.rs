//! A card's Newtonian liquid in a laminar circular pipe.
//!
//! [`MaterialLaminarPipe`] is the Hagen-Poiseuille resistance of fully
//! developed laminar flow, `Q = G (p_a - p_b)` with `G = pi r^4 / (8 mu L)`.
//! Its dynamic viscosity `mu` is a material card's, at the temperature it is
//! asked at: the card's dynamic-viscosity claim, or its kinematic-viscosity
//! claim times its density at the same point. All of the viscous power
//! `G (p_a - p_b)^2` is heat. Every evaluation resolves the SAME pinned claims
//! through the card again ([`ResolvedMaterialStatePoint::at_query_point`]),
//! so a temperature outside the declared span or the card's own support is
//! refused, never extrapolated. The temperature slope is that of the card's
//! pinned line on the segment the temperature lies in: the right segment at
//! an interior knot and the left one at the last knot, as the Ohmic node
//! takes its resistance slope. With a kinematic claim the slope is the product
//! rule's, `d(rho nu)/dT = rho' nu + rho nu'`, on the same segments.
//!
//! [`LaminarPipeComponent`] (feature `component`) puts the pipe behind two
//! hydraulic ports and one heat port, the ports of fs-material's Carreau pipe.
//!
//! This is a quasisteady, incompressible resistance with one temperature for
//! the whole liquid: no fluid storage, inertia, entrance length, transition to
//! turbulence (laminar flow is the caller's assumption; no Reynolds number is
//! checked), wall slip, pressure dependence of the viscosity, thermal expansion
//! of the pipe or heat transport along it. The card's uncertainty is not
//! carried.

use std::sync::Arc;

use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_qty::{Density, Temperature};

use crate::liquid::{DYNAMIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_DIMS};
use crate::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ResolvedScalarProperty,
    ScalarAdmissibility, is_fixed_external_coordinate, pinned_temperature_line,
};

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

/// The slope of a pinned `(T, value)` line on the segment `temperature_k`
/// lies in: the right segment at an interior knot, the left one at the last.
fn segment_slope(line: &[(f64, f64)], temperature_k: f64) -> f64 {
    let segment = line
        .partition_point(|&(t, _)| t <= temperature_k)
        .saturating_sub(1)
        .min(line.len() - 2);
    let [left, right] = [line[segment], line[segment + 1]];
    (right.1 - left.1) / (right.0 - left.0)
}

/// The viscosity and the density, if any, of a pipe's source. Refused: a
/// property that is neither, two of either, no viscosity, a dynamic viscosity
/// with a density, and a kinematic one without.
fn classify(
    original: &ResolvedMaterialStatePoint,
) -> Result<(&ResolvedScalarProperty, Option<&ResolvedScalarProperty>), MaterialStatePointError> {
    let (mut viscosity, mut density) = (None, None);
    for property in original.properties() {
        let dims = property.requirement().dims();
        let slot = if dims == DYNAMIC_VISCOSITY_DIMS || dims == KINEMATIC_VISCOSITY_DIMS {
            &mut viscosity
        } else if dims == Density::DIMS {
            &mut density
        } else {
            return Err(invalid(
                "laminar pipe source holds a property that is neither a viscosity nor a density",
            ));
        };
        if slot.replace(property).is_some() {
            return Err(invalid(
                "laminar pipe source holds two viscosities or two densities",
            ));
        }
    }
    let viscosity = viscosity.ok_or_else(|| invalid("laminar pipe source holds no viscosity"))?;
    match (
        viscosity.requirement().dims() == DYNAMIC_VISCOSITY_DIMS,
        density.is_some(),
    ) {
        (true, true) => Err(invalid(
            "laminar pipe with a dynamic viscosity takes no density",
        )),
        (false, false) => Err(invalid(
            "laminar pipe with a kinematic viscosity needs the card's density",
        )),
        _ => Ok((viscosity, density)),
    }
}

/// The card's selected claim for `property`, pinned along `axis` from the
/// source coordinates at the range's two ends, its knots checked against the
/// range, the original admissibility and a representable slope.
fn pinned_line(
    card: &MaterialCard,
    property: &ResolvedScalarProperty,
    axis: &str,
    [from, to]: [&QueryPoint; 2],
    [low, high]: [f64; 2],
) -> Result<(PropertyKey, Vec<(f64, f64)>), MaterialStatePointError> {
    if !property.fixed_external_support(Some(axis), false, true) {
        return Err(invalid(
            "laminar pipe requires one thermal coordinate and fixed external source conditions",
        ));
    }
    let selected = property.answer().receipt.selected;
    let claim = card
        .claims()
        .claim(selected)
        .ok_or_else(|| invalid("laminar pipe source claim is absent"))?;
    let line = pinned_temperature_line(card, &claim.key, axis, from, to, selected)?;
    if line.len() < 2
        || line.first().map(|p| p.0) != Some(low)
        || line.last().map(|p| p.0) != Some(high)
    {
        return Err(invalid(
            "laminar pipe source line does not span the complete range",
        ));
    }
    let admissibility = property.requirement().admissibility();
    for &(_, value) in &line {
        let admitted = value.is_finite()
            && match admissibility {
                ScalarAdmissibility::Finite => true,
                ScalarAdmissibility::StrictlyPositive => value > 0.0,
                ScalarAdmissibility::NonNegative => value >= 0.0,
                ScalarAdmissibility::OpenInterval { lower, upper } => {
                    lower < value && value < upper
                }
            };
        if !admitted {
            return Err(MaterialStatePointError::OutsideConsumerDomain {
                property: claim.key.name().to_owned(),
                value,
                admissibility,
            });
        }
    }
    for pair in line.windows(2) {
        let width = pair[1].0 - pair[0].0;
        let slope = (pair[1].1 - pair[0].1) / width;
        if !width.is_finite() || width <= 0.0 || !slope.is_finite() {
            return Err(invalid("laminar pipe source slope is unrepresentable"));
        }
    }
    Ok((claim.key.clone(), line))
}

/// How the card states the liquid's viscosity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ViscosityBasis {
    /// A dynamic-viscosity claim [Pa s].
    Dynamic,
    /// A kinematic-viscosity claim [m²/s], times the card's density there.
    Kinematic,
}

/// A laminar circular pipe of fixed radius and length whose liquid's
/// viscosity is a material card's at the temperature it is asked at.
#[derive(Clone, Debug)]
pub struct MaterialLaminarPipe {
    card: Arc<MaterialCard>,
    original: ResolvedMaterialStatePoint,
    viscosity: (PropertyKey, Vec<(f64, f64)>),
    density: Option<(PropertyKey, Vec<(f64, f64)>)>,
    temperature_axis: String,
    reference_temperature_k: f64,
    temperature_range_k: [f64; 2],
    radius_m: f64,
    length_m: f64,
    /// `pi r^4 / (8 L)` [m³].
    geometry_m3: f64,
}

impl MaterialLaminarPipe {
    /// Prepare a pipe over `temperature_range_k` from `original`, a state
    /// point of `card` that holds a dynamic viscosity alone, or a kinematic
    /// viscosity and a density, resolved at a point whose `temperature_axis`
    /// is an absolute temperature inside the range.
    ///
    /// Every other coordinate of the point must be a fixed external condition
    /// of the claims (a dimensionless provenance flag pinned to one value is
    /// held, as the thermal owner holds it). Each claim's pinned line over the
    /// range is read once for the slopes, and the card is resolved at every
    /// knot of those lines, so a claim the range leaves, or a viscosity that is
    /// not positive somewhere in it, refuses here and not in a later step.
    ///
    /// # Errors
    /// Another card, a property that is neither a viscosity nor a density, a
    /// dynamic viscosity with a density or a kinematic one without, a
    /// temperature coordinate that is not an absolute temperature, a range
    /// that is not positive or does not hold the original temperature, a
    /// coordinate that is not a fixed external condition, a geometry that is
    /// not positive and finite, or the card's refusal at a knot.
    pub fn new(
        card: Arc<MaterialCard>,
        original: ResolvedMaterialStatePoint,
        temperature_axis: &str,
        temperature_range_k: [f64; 2],
        radius_m: f64,
        length_m: f64,
    ) -> Result<Self, MaterialStatePointError> {
        let geometry_m3 = core::f64::consts::PI * radius_m.powi(4) / (8.0 * length_m);
        if !(radius_m.is_finite()
            && radius_m > 0.0
            && length_m.is_finite()
            && length_m > 0.0
            && geometry_m3.is_finite()
            && geometry_m3 > 0.0)
        {
            return Err(invalid(
                "laminar pipe radius and length must be positive and finite",
            ));
        }
        let (viscosity, density) = classify(&original)?;
        let context = original.query();
        let reference_temperature_k = *context
            .axes()
            .get(temperature_axis)
            .ok_or_else(|| invalid("laminar pipe source temperature is absent"))?;
        let quantity = context.axis_quantities().get(temperature_axis).copied();
        if quantity.is_some_and(|q| q.dims() != Temperature::DIMS)
            || !is_fixed_external_coordinate(
                temperature_axis,
                quantity,
                Some(reference_temperature_k),
            )
        {
            return Err(invalid(
                "laminar pipe needs an original absolute-temperature coordinate",
            ));
        }
        let [low, high] = temperature_range_k;
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || low >= high
            || !(low..=high).contains(&reference_temperature_k)
        {
            return Err(invalid(
                "laminar pipe range must contain the original positive temperature",
            ));
        }
        let lower = context
            .clone()
            .with(temperature_axis, low)
            .map_err(|_| invalid("laminar pipe lower source coordinate is invalid"))?;
        let upper = context
            .with(temperature_axis, high)
            .map_err(|_| invalid("laminar pipe upper source coordinate is invalid"))?;
        let line = |property| {
            pinned_line(
                &card,
                property,
                temperature_axis,
                [&lower, &upper],
                temperature_range_k,
            )
        };
        let viscosity = line(viscosity)?;
        let density = density.map(line).transpose()?;
        let pipe = Self {
            card,
            original,
            viscosity,
            density,
            temperature_axis: temperature_axis.to_owned(),
            reference_temperature_k,
            temperature_range_k,
            radius_m,
            length_m,
            geometry_m3,
        };
        // The card itself at every knot of every line: its own support, the
        // original admissibility, and a positive, representable viscosity.
        let knots = pipe
            .viscosity
            .1
            .iter()
            .chain(pipe.density.iter().flat_map(|(_, line)| line))
            .map(|&(temperature, _)| temperature);
        for temperature in knots {
            pipe.conductance_and_slope(temperature)?;
        }
        Ok(pipe)
    }

    /// The liquid's dynamic viscosity [Pa s] at `temperature_k` and its
    /// temperature slope [Pa s/K], both from the card.
    ///
    /// # Errors
    /// A temperature outside the declared range, the card's refusal there, or
    /// a viscosity that is not positive and representable.
    pub fn viscosity_and_slope(
        &self,
        temperature_k: f64,
    ) -> Result<(f64, f64), MaterialStatePointError> {
        let [low, high] = self.temperature_range_k;
        if !(temperature_k.is_finite() && low <= temperature_k && temperature_k <= high) {
            return Err(invalid(
                "laminar pipe temperature outside its declared range",
            ));
        }
        let query = self
            .original
            .query()
            .with(&self.temperature_axis, temperature_k)
            .map_err(|_| invalid("laminar pipe temperature coordinate is invalid"))?;
        let state = self.original.at_query_point(&self.card, &query)?;
        let value = |key: &PropertyKey| {
            state
                .property_by_key(key)
                .map(ResolvedScalarProperty::value_si)
                .ok_or_else(|| invalid("laminar pipe refreshed source lost a property"))
        };
        let (key, line) = &self.viscosity;
        let stated = value(key)?;
        let stated_slope = segment_slope(line, temperature_k);
        let (viscosity, slope) = match &self.density {
            None => (stated, stated_slope),
            Some((key, line)) => {
                let density = value(key)?;
                let density_slope = segment_slope(line, temperature_k);
                (
                    density * stated,
                    density_slope.mul_add(stated, density * stated_slope),
                )
            }
        };
        if !(viscosity.is_finite() && viscosity > 0.0 && slope.is_finite()) {
            return Err(invalid(
                "laminar pipe viscosity or its slope is unrepresentable",
            ));
        }
        Ok((viscosity, slope))
    }

    /// The pipe's hydraulic conductance `G = pi r^4 / (8 mu L)` [m³/(s Pa)]
    /// at `temperature_k` and its temperature slope `-G mu' / mu`.
    ///
    /// # Errors
    /// As [`Self::viscosity_and_slope`], or an unrepresentable conductance.
    pub fn conductance_and_slope(
        &self,
        temperature_k: f64,
    ) -> Result<(f64, f64), MaterialStatePointError> {
        let (viscosity, slope) = self.viscosity_and_slope(temperature_k)?;
        let conductance = self.geometry_m3 / viscosity;
        let conductance_slope = -conductance * slope / viscosity;
        if !(conductance.is_finite() && conductance > 0.0 && conductance_slope.is_finite()) {
            return Err(invalid("laminar pipe conductance is unrepresentable"));
        }
        Ok((conductance, conductance_slope))
    }

    /// How the card states the viscosity.
    #[must_use]
    pub const fn basis(&self) -> ViscosityBasis {
        if self.density.is_some() {
            ViscosityBasis::Kinematic
        } else {
            ViscosityBasis::Dynamic
        }
    }

    /// The card the pipe resolves.
    #[must_use]
    pub fn card(&self) -> &MaterialCard {
        &self.card
    }

    /// The state point the pipe was prepared from.
    #[must_use]
    pub const fn original(&self) -> &ResolvedMaterialStatePoint {
        &self.original
    }

    /// The coordinate the heat port's temperature sets.
    #[must_use]
    pub fn temperature_axis(&self) -> &str {
        &self.temperature_axis
    }

    /// The original state point's temperature [K].
    #[must_use]
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }

    /// The declared range the pipe answers in [K].
    #[must_use]
    pub const fn temperature_range_k(&self) -> [f64; 2] {
        self.temperature_range_k
    }

    /// The bore radius [m].
    #[must_use]
    pub const fn radius_m(&self) -> f64 {
        self.radius_m
    }

    /// The length [m].
    #[must_use]
    pub const fn length_m(&self) -> f64 {
        self.length_m
    }
}

#[cfg(feature = "component")]
pub use component::LaminarPipeComponent;

#[cfg(feature = "component")]
mod component {
    use core::num::NonZeroUsize;

    use fs_exec::Cx;
    use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
    use fs_ports::component::{
        Class, Component, EndOfStep, PortRole, PortSpec, PortTangent, Proposal, Receipt,
    };
    use fs_ports::{PortKind, StableId};

    use super::{MaterialLaminarPipe, ViscosityBasis, invalid};
    use crate::state_point::MaterialStatePointError;

    fn poll(cx: &Cx<'_>) -> Result<(), MaterialStatePointError> {
        cx.checkpoint()
            .map_err(|_| invalid("laminar pipe cancelled"))
    }

    fn values(values: &[f64]) -> Result<[f64; 3], MaterialStatePointError> {
        match values {
            [a, b, c] if values.iter().all(|value| value.is_finite()) => Ok([*a, *b, *c]),
            _ => Err(invalid("laminar pipe needs three finite port values")),
        }
    }

    /// A card's laminar pipe behind two hydraulic ports and one heat port.
    ///
    /// Inputs, in port order, are the pressures at `a` and `b` [Pa] and the
    /// liquid's temperature at `heat` [K]. Outputs are the volume flows `Q`
    /// and `-Q` into the pipe [m³/s] and `-G (p_a - p_b)^2` [W]: all of the
    /// viscous power leaves through the heat port. Each proposal resolves the
    /// card again at the heat port's temperature. The state is `()`.
    ///
    /// Refusals are the pipe's own [`MaterialStatePointError`]: the card's
    /// support, the declared range, and `InvalidDerived` naming a port value,
    /// timestep or slot count the pipe cannot take, or a cancelled caller. No
    /// output is written on a refusal.
    #[derive(Clone, Debug)]
    pub struct LaminarPipeComponent {
        pipe: MaterialLaminarPipe,
        ports: [PortSpec; 3],
    }

    impl LaminarPipeComponent {
        /// Put `pipe` behind its ports.
        #[must_use]
        pub fn new(pipe: MaterialLaminarPipe) -> Self {
            let hydraulic = |name| {
                PortSpec::power(
                    StableId::new(name).expect("static port name"),
                    PortKind::FluidPressureFlux,
                    PortRole::EffortIn,
                    NonZeroUsize::MIN,
                )
                .expect("scalar hydraulic power port")
            };
            Self {
                pipe,
                ports: [
                    hydraulic("a"),
                    hydraulic("b"),
                    PortSpec::heat(
                        StableId::new("heat").expect("static port name"),
                        PortRole::EffortIn,
                    )
                    .expect("thermal outlet"),
                ],
            }
        }

        /// The pipe behind the ports.
        #[must_use]
        pub const fn pipe(&self) -> &MaterialLaminarPipe {
            &self.pipe
        }

        /// The Jacobian of `[Q, -Q, -P]` in `[p_a, p_b, T]`.
        fn jacobian(&self, inputs: &[f64]) -> Result<[[f64; 3]; 3], MaterialStatePointError> {
            let [a, b, kelvin] = values(inputs)?;
            let (conductance, slope) = self.pipe.conductance_and_slope(kelvin)?;
            let drop = a - b;
            let flow_temperature = slope * drop;
            let heat_pressure = 2.0 * conductance * drop;
            let heat_temperature = flow_temperature * drop;
            let matrix = [
                [conductance, -conductance, flow_temperature],
                [-conductance, conductance, -flow_temperature],
                [-heat_pressure, heat_pressure, -heat_temperature],
            ];
            if matrix.iter().flatten().any(|value| !value.is_finite()) {
                return Err(invalid("laminar pipe tangent is unrepresentable"));
            }
            Ok(matrix)
        }
    }

    impl Component for LaminarPipeComponent {
        type State = ();
        type Error = MaterialStatePointError;

        fn ports(&self) -> &[PortSpec] {
            &self.ports
        }

        fn class(&self) -> Class {
            // Smooth inside a segment of the card's lines.
            Class::Smooth
        }

        fn stored_energy_j(&self, _state: &()) -> Result<f64, MaterialStatePointError> {
            Ok(0.0)
        }

        const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
            owner: "fs-material::fluid_transport::laminar_pipe::LaminarPipeComponent",
            state_family: "stateless-card-laminar-pipe",
            schema_grammar: "unit state on the card's pinned viscosity claims, range and geometry",
            codec_grammar: "zero u64 words",
            codec_version: 1,
        });

        fn encode_state(&self, _state: &(), _words: &mut Vec<u64>) -> bool {
            true
        }

        fn decode_state(&self, words: &[u64]) -> Option<()> {
            words.is_empty().then_some(())
        }

        fn identity_words(&self, words: &mut Vec<u64>) {
            // The card, the selected claims and the point they were resolved
            // at; the range, the reference temperature and the geometry. The
            // temperature of each step is not bound: every proposal resolves
            // the card again at the heat port's.
            let hash = |bytes: &[u8; 32]| {
                bytes
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|chunk| u64::from_le_bytes(*chunk))
                    .collect::<Vec<_>>()
            };
            words.extend(hash(self.pipe.card().content_hash().as_bytes()));
            for property in self.pipe.original().properties() {
                words.extend(hash(property.answer().receipt.content_hash().as_bytes()));
            }
            let axis = self.pipe.temperature_axis().as_bytes();
            words.push(axis.len() as u64);
            words.extend(axis.iter().map(|&byte| u64::from(byte)));
            words.push(self.pipe.reference_temperature_k().to_bits());
            words.extend(self.pipe.temperature_range_k().map(f64::to_bits));
            words.extend([
                self.pipe.radius_m().to_bits(),
                self.pipe.length_m().to_bits(),
                match self.pipe.basis() {
                    ViscosityBasis::Dynamic => 0,
                    ViscosityBasis::Kinematic => 1,
                },
            ]);
        }

        fn end_of_step(&self, _port: usize) -> EndOfStep {
            // An algebraic, held-input response: no flow history from which an
            // endpoint relation could be reconstructed.
            EndOfStep::Untold
        }

        fn propose(
            &self,
            _old: &(),
            dt_s: f64,
            inputs: &[f64],
            new: &mut (),
            outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<Proposal, MaterialStatePointError> {
            poll(cx)?;
            if outputs.len() != 3 {
                return Err(invalid("laminar pipe needs three output slots"));
            }
            if !(dt_s.is_finite() && dt_s > 0.0) {
                return Err(invalid("laminar pipe timestep must be positive and finite"));
            }
            let [a, b, kelvin] = values(inputs)?;
            let (conductance, _) = self.pipe.conductance_and_slope(kelvin)?;
            let drop = a - b;
            let flow = conductance * drop;
            let power = flow * drop;
            let dissipated_j = dt_s * power;
            let energy_tolerance_j = 64.0
                * f64::EPSILON
                * (dt_s * (flow * a).abs() + dt_s * (flow * b).abs() + dissipated_j);
            let entropy = dissipated_j / kelvin;
            if ![flow, power, dissipated_j, energy_tolerance_j, entropy]
                .into_iter()
                .all(f64::is_finite)
            {
                return Err(invalid("laminar pipe energy account is unrepresentable"));
            }
            poll(cx)?;
            *new = ();
            outputs.copy_from_slice(&[flow, -flow, -power]);
            Ok(Proposal::Step(Receipt {
                stored_change_j: 0.0,
                dissipated_j,
                numerical_loss_j: 0.0,
                entropy_generated_j_per_k: entropy,
                energy_tolerance_j,
                residual: 0.0,
                tolerance: 0.0,
            }))
        }
    }

    impl PortTangent for LaminarPipeComponent {
        fn jvp(
            &self,
            _old: &(),
            _dt_s: f64,
            inputs: &[f64],
            d_inputs: &[f64],
            d_outputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), MaterialStatePointError> {
            poll(cx)?;
            if d_outputs.len() != 3 {
                return Err(invalid("laminar pipe needs three tangent output slots"));
            }
            let direction = values(d_inputs)?;
            let matrix = self.jacobian(inputs)?;
            let mut next = [0.0; 3];
            for (value, row) in next.iter_mut().zip(matrix) {
                *value = row.iter().zip(direction).map(|(j, d)| j * d).sum();
            }
            d_outputs.copy_from_slice(&next);
            Ok(())
        }

        fn vjp(
            &self,
            _old: &(),
            _dt_s: f64,
            inputs: &[f64],
            bar_outputs: &[f64],
            bar_inputs: &mut [f64],
            cx: &Cx<'_>,
        ) -> Result<(), MaterialStatePointError> {
            poll(cx)?;
            if bar_inputs.len() != 3 {
                return Err(invalid("laminar pipe needs three tangent input slots"));
            }
            let direction = values(bar_outputs)?;
            let matrix = self.jacobian(inputs)?;
            let mut next = [0.0; 3];
            for (column, value) in next.iter_mut().enumerate() {
                *value = matrix
                    .iter()
                    .zip(direction)
                    .map(|(row, d)| row[column] * d)
                    .sum();
            }
            bar_inputs.copy_from_slice(&next);
            Ok(())
        }
    }
}
