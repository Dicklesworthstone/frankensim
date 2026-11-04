//! Coexisting liquid phases enter the original open two-fluid Washburn owner.
//!
//! Surface tension is explicitly authored or derived from the original
//! source-bound square-gradient interface. Contact angle remains caller-authored
//! for the ACTUAL coexisting pair and substrate at this temperature. Neither
//! path borrows a liquid/gas surface or pure-liquid wetting card. An equilibrium
//! with one miscible phase has no liquid-liquid meniscus and is refused before
//! any flow or work is evaluated.

use super::{MaterialCapillaryError, MaterialMeniscusInput, signed_mass};
use crate::capillary::{
    CircularMeniscusInput, TwoFluidCircularMeniscusInput, TwoFluidCircularMeniscusStep,
    step_newtonian_two_fluid_circular_meniscus,
};
use fs_material::liquid::solution::gradient::interface::ResolvedBinaryLiquidPlanarInterface;
use fs_material::liquid::solution::{ResolvedBinaryLiquidEquilibrium, ResolvedBinaryLiquidPhase};

/// Accepted open-reservoir displacement, retaining the original source phases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidMeniscusStep<'a> {
    /// Original equilibrium and source receipts; its finite inventory is unchanged.
    pub equilibrium: &'a ResolvedBinaryLiquidEquilibrium,
    /// Invading phase in the equilibrium owner's composition order (0 or 1).
    pub invading_phase: usize,
    /// Geometry, current position, applied pressure and requested duration.
    pub input: MaterialMeniscusInput,
    /// Explicit pair tension [N/m]; zero explicitly neglects capillary work.
    pub surface_tension_n_m: f64,
    /// Explicit equilibrium wall angle, measured through the invading phase.
    pub contact_angle_cosine: f64,
    /// Existing exact end-event, signed volume, work and two-fluid loss account.
    pub discrete: TwoFluidCircularMeniscusStep,
    /// Signed component amounts entering at the inlet [mol], in species order.
    pub invading_component_inflow_mol: [f64; 2],
    /// Signed component amounts leaving at the outlet [mol], in species order.
    pub displaced_component_outflow_mol: [f64; 2],
    /// Net tube inventory change: inlet amount minus outlet amount [mol].
    pub tube_component_change_mol: [f64; 2],
    /// Signed invading mass entering at the inlet [kg].
    pub invading_mass_inflow_kg: f64,
    /// Signed displaced mass leaving at the outlet [kg].
    pub displaced_mass_outflow_kg: f64,
}

/// Sharp-meniscus approximation retaining the SAME planar-interface owner and
/// the complete original two-fluid flow, work and signed component account.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BinaryLiquidPlanarMeniscusStep<'a> {
    /// Source phases, original bulk law, kappa and interface quadrature receipts.
    pub planar_interface: &'a ResolvedBinaryLiquidPlanarInterface,
    /// Characteristic central-slope interface width divided by tube radius.
    pub interface_width_radius_ratio: f64,
    /// Explicit caller allowance used to qualify the sharp-meniscus reduction.
    pub maximum_interface_width_radius_ratio: f64,
    /// Existing accepted step with sigma supplied by this planar interface.
    pub meniscus: BinaryLiquidMeniscusStep<'a>,
}

/// Advance the same coexisting liquids using their derived planar tension.
///
/// The caller supplies the wall angle and the maximum allowed characteristic
/// interface width/radius, a finite number in (0,1]. The central-slope width is
/// the composition gap divided by the maximum planar slope; it is not the full
/// extent of the profile's infinite tails. Selecting a suitably small allowance
/// is part of qualifying this planar-to-curved, sharp-meniscus approximation.
/// The ratio check itself does not bound curvature or contact-line errors.
///
/// No wall angle is inferred from bulk miscibility, and the calculated tension
/// is not promoted to a measured material claim. All open-reservoir, frozen-
/// coefficient and pressure-independent-coexistence assumptions of
/// [`step_binary_liquid_horizontal`] still apply. This adapter allocates no memory.
///
/// # Errors
/// Refuses an invalid width allowance, unrepresentable geometry ratio or an
/// interface too broad for the caller's allowance before evaluating any flow.
/// The existing two-fluid owner then admits the original phase and flow inputs.
pub fn step_binary_liquid_planar_horizontal(
    planar_interface: &ResolvedBinaryLiquidPlanarInterface,
    invading_phase: usize,
    contact_angle_cosine: f64,
    maximum_interface_width_radius_ratio: f64,
    input: MaterialMeniscusInput,
) -> Result<BinaryLiquidPlanarMeniscusStep<'_>, MaterialCapillaryError> {
    if !maximum_interface_width_radius_ratio.is_finite()
        || maximum_interface_width_radius_ratio <= 0.0
        || maximum_interface_width_radius_ratio > 1.0
        || !input.radius_m.is_finite()
        || input.radius_m <= 0.0
    {
        return Err(MaterialCapillaryError::Unsupported(
            "positive tube radius and finite interface-width/radius allowance in (0,1] required",
        ));
    }
    let interface = planar_interface.interface();
    let ratio = interface.central_slope_width_m() / input.radius_m;
    if !ratio.is_finite() || ratio <= 0.0 || ratio > maximum_interface_width_radius_ratio {
        return Err(MaterialCapillaryError::Unsupported(
            "planar interface width exceeds the admitted sharp-meniscus scale or is unrepresentable",
        ));
    }
    let meniscus = step_binary_liquid_horizontal(
        planar_interface.equilibrium(),
        invading_phase,
        interface.surface_tension_n_m(),
        contact_angle_cosine,
        input,
    )?;
    Ok(BinaryLiquidPlanarMeniscusStep {
        planar_interface,
        interface_width_radius_ratio: ratio,
        maximum_interface_width_radius_ratio,
        meniscus,
    })
}

/// Advance a horizontal meniscus between the two actual equilibrium liquids.
///
/// Both phases are homogeneous, Newtonian and incompressible with frozen source
/// coefficients. Their reservoirs maintain the supplied compositions and inlet/
/// outlet pressures; the finite equilibrium inventory is used to identify the
/// phases and concentrations, not as a reservoir that is silently depleted.
/// Density and viscosity sources must admit static fixed external conditions;
/// rate/frequency dependence or applicability does not become Newtonian merely
/// because a source was evaluated once. Held external axes must be typed.
/// Reverse motion returns negative invading inflow and negative displaced
/// outflow. These signed transfers close the tube's component and mass balances.
///
/// The declared bulk coexistence is used despite the capillary pressure jump;
/// pressure-dependent coexistence, interfacial mass exchange, diffusion,
/// nucleation, films, inertia and dynamic contact angles are outside this law.
/// The caller must supply pair/substrate interface data valid under these
/// assumptions. This bounded scalar adapter allocates no memory.
///
/// # Errors
/// Refuses one-phase equilibrium, an invalid phase index, invalid authored
/// interface/flow inputs and unrepresentable transfers, without mutating state.
pub fn step_binary_liquid_horizontal(
    equilibrium: &ResolvedBinaryLiquidEquilibrium,
    invading_phase: usize,
    surface_tension_n_m: f64,
    contact_angle_cosine: f64,
    input: MaterialMeniscusInput,
) -> Result<BinaryLiquidMeniscusStep<'_>, MaterialCapillaryError> {
    let [first, second] = equilibrium.phases() else {
        return Err(MaterialCapillaryError::Unsupported(
            "a liquid-liquid meniscus requires two coexisting phases; this solution is miscible",
        ));
    };
    let (invading, displaced) = match invading_phase {
        0 => (first, second),
        1 => (second, first),
        _ => {
            return Err(MaterialCapillaryError::Unsupported(
                "binary invading phase index must be 0 or 1",
            ));
        }
    };
    for phase in [invading, displaced] {
        let liquid = phase.liquid();
        let requirements = liquid.requirements().properties();
        for requirement in [&requirements[0], &requirements[3]] {
            if !liquid
                .material()
                .property(requirement.name())
                .is_some_and(|property| property.is_static_at_fixed_external_conditions())
            {
                return Err(MaterialCapillaryError::Unsupported(
                    "binary capillary density and viscosity require static fixed external conditions",
                ));
            }
        }
    }
    let discrete = step_newtonian_two_fluid_circular_meniscus(TwoFluidCircularMeniscusInput {
        meniscus: CircularMeniscusInput {
            dynamic_viscosity_pa_s: invading.liquid().dynamic_viscosity_pa_s(),
            radius_m: input.radius_m,
            tube_length_m: input.tube_length_m,
            wet_length_m: input.wet_length_m,
            surface_tension_n_m,
            contact_angle_cosine,
            applied_pressure_pa: input.applied_pressure_pa,
            duration_s: input.duration_s,
        },
        displaced_dynamic_viscosity_pa_s: displaced.liquid().dynamic_viscosity_pa_s(),
    })?;
    let volume = discrete.meniscus.transferred_volume_m3;
    let incoming = component_transfer(invading, volume)?;
    let outgoing = component_transfer(displaced, volume)?;
    let net = [incoming[0] - outgoing[0], incoming[1] - outgoing[1]];
    if net.iter().any(|n| !n.is_finite()) {
        return Err(MaterialCapillaryError::NonFiniteComponentTransfer);
    }
    Ok(BinaryLiquidMeniscusStep {
        equilibrium,
        invading_phase,
        input,
        surface_tension_n_m,
        contact_angle_cosine,
        discrete,
        invading_component_inflow_mol: incoming,
        displaced_component_outflow_mol: outgoing,
        tube_component_change_mol: net,
        invading_mass_inflow_kg: signed_mass(invading.liquid().density_kg_m3(), volume)?,
        displaced_mass_outflow_kg: signed_mass(displaced.liquid().density_kg_m3(), volume)?,
    })
}

fn component_transfer(
    phase: &ResolvedBinaryLiquidPhase,
    volume_m3: f64,
) -> Result<[f64; 2], MaterialCapillaryError> {
    if volume_m3 == 0.0 {
        return Ok([0.0; 2]);
    }
    let amounts = phase.component_moles();
    let transfer = amounts.map(|amount| (amount / phase.volume_m3()) * volume_m3);
    for (amount, moved) in amounts.into_iter().zip(transfer) {
        if !moved.is_finite() || (amount > 0.0 && moved == 0.0) {
            return Err(MaterialCapillaryError::NonFiniteComponentTransfer);
        }
    }
    Ok(transfer)
}
