//! A source convection law on the original material/phase/finite-gas timeline.
//!
//! The existing bulk conductivity and phase correction are held over the half-cell
//! at each epoch. The common material surface solver closes this local path
//! against actual gas film properties. This does not resolve a phase front or
//! a conductivity profile inside that half-cell, nor add thermal storage.

use super::*;
use fs_evidence::Evidence;
use fs_material::thermal::conduction::{
    ConductivityCurve, ConvectiveLayerError, LayerConductionConfig, ReservoirLayerGeometry,
};

/// Gas-side h [W/(m² K)], independent of the higher correlation crate.
pub type FiniteGasHeatTransferCoefficient = fs_qty::Qty<0, 1, -3, -1, 0, 0>;

/// Stationary physical patch with an explicitly authored material-to-gas normal.
/// The original exposed-cell coverage is still checked by the native owner.
/// Scalar contact area alone cannot establish orientation; callers must supply
/// the physical unit normal and retain the same geometry throughout the split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteGasConvectionContact {
    /// Positive area and material distance; gas thickness must be ZERO.
    pub contact: FiniteGasContact,
    pub outward_normal: [f64; 3],
}

/// Actual surface trial and the original gas at its current EOS pressure.
#[derive(Clone, Copy, Debug)]
pub struct FiniteGasConvectionQuery<'a> {
    pub contact: FiniteGasConvectionContact,
    pub condensed_temperature_k: f64,
    pub surface_temperature_k: f64,
    pub bulk: &'a ResolvedThermochemicalGasState,
    /// Same source claims at (surface T + bulk T)/2, never a second gas store.
    pub film: &'a ResolvedThermochemicalGasState,
}

/// Repeatable constitutive source, with geometry/regime admission at EVERY trial.
/// `None` denotes only exact surface/bulk equilibrium. Driven h must be finite
/// and positive. Rejected trials and transactional clones may repeat queries;
/// providers must not advance physical history, integrate heat or replace a
/// refused source. All physical lengths and relative flow remain provider data.
pub trait FiniteGasConvection: core::fmt::Debug {
    fn coefficient(
        &self,
        query: FiniteGasConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FiniteGasHeatTransferCoefficient>>, MaterialError>;
}

/// Accepted nonlinear surface and original gas-side correlation evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasConvectionSample {
    pub outward_normal: [f64; 3],
    pub surface_temperature_k: f64,
    /// Positive from the original material half-cell toward its gas surface.
    pub condensed_power_to_gas_w: f64,
    /// Material integral-coordinate heat minus actual gas heat [W].
    pub surface_heat_residual_w: f64,
    pub iterations: usize,
    pub coefficient: Option<Evidence<FiniteGasHeatTransferCoefficient>>,
}

#[derive(Clone, Copy)]
pub(super) struct ConvectionBinding<'a> {
    contacts: &'a [FiniteGasConvectionContact],
    provider: &'a dyn FiniteGasConvection,
    config: LayerConductionConfig,
}

impl ConvectionBinding<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn sample<E>(
        self,
        index: usize,
        gas: &IsochoricGasParcel,
        specific_enthalpy_j_kg: f64,
        cell_k: f64,
        conductivity: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FiniteGasContactSample, FiniteGasFlowError<E>> {
        check()?;
        let contact = *self.contacts.get(index)
            .ok_or(invalid("finite gas convection lost its original contact"))?;
        let bulk = gas.state();
        let bulk_k = bulk.state().temperature;
        let area = contact.contact.area_m2;
        let distance = contact.contact.distances_m[0];
        if !conductivity.is_finite() || conductivity < 0.0 {
            return Err(DepositionError::Arithmetic.into());
        }
        let mut film_at = |surface: f64| -> Result<_, FiniteGasFlowError<E>> {
            let point = gas.point().clone()
                .with("temperature", 0.5 * surface + 0.5 * bulk_k)
                .map_err(|_| invalid("invalid finite gas convection film temperature"))?;
            let film = gas.model().evaluate(&point).map_err(GasEnergyError::from)?;
            if !same_sources(bulk, &film) {
                return Err(invalid("finite gas surface left its original source claims").into());
            }
            let coefficient = self.provider.coefficient(FiniteGasConvectionQuery {
                contact,
                condensed_temperature_k: cell_k,
                surface_temperature_k: surface,
                bulk,
                film: &film,
            }).map_err(|source| MaterialFlowError::Material {
                cell: contact.contact.cell, source,
            })?;
            let conductance = match coefficient.as_ref() {
                Some(h) if h.value.value().is_finite() && h.value.value() > 0.0 => {
                    let g = finite(area * h.value.value())?;
                    if g <= 0.0 { return Err(DepositionError::Arithmetic.into()); }
                    g
                }
                None if surface == bulk_k => 0.0,
                _ => return Err(invalid("finite gas convection requires positive finite h or exact equilibrium").into()),
            };
            let power = finite(conductance * (surface - bulk_k))?;
            Ok((power, conductance, (film, coefficient)))
        };

        let (surface, solid_power, power, conductance, residual, iterations, film, coefficient) =
            if conductivity == 0.0 || cell_k == bulk_k {
                // A declared insulator has zero heat and its exterior surface
                // equilibrates with gas. Exact equilibrium never needs a fake
                // finite-width conductivity chart or an invented finite-Ra h.
                let (_, gas_g, (film, coefficient)) = film_at(bulk_k)?;
                let g = if conductivity == 0.0 || gas_g == 0.0 {
                    0.0
                } else {
                    let resistance = finite(distance / conductivity + area / gas_g)?;
                    let g = finite(area / resistance)?;
                    if resistance <= 0.0 || g <= 0.0 {
                        return Err(DepositionError::Arithmetic.into());
                    }
                    g
                };
                (bulk_k, 0.0, 0.0, g, 0.0, 0, film, coefficient)
            } else {
                // This is the native epoch's HELD scalar k, including its
                // existing phase treatment. These two knots declare that
                // finite-volume approximation; they are not fabricated source
                // samples or an extrapolated pure-phase material chart.
                let curve = ConductivityCurve::new(vec![
                    [cell_k.min(bulk_k), conductivity],
                    [cell_k.max(bulk_k), conductivity],
                ])?;
                let response = ReservoirLayerGeometry::planar(area, distance)?
                    .evaluate_convective_material_to_bath(
                        &curve, [cell_k, bulk_k], self.config, &mut film_at, &mut *check,
                    ).map_err(|error| match error {
                        ConvectiveLayerError::Thermal(error) => error.into(),
                        ConvectiveLayerError::Film(error) => error,
                        ConvectiveLayerError::SurfaceBalance { residual_w, iterations } => {
                            FiniteGasFlowError::SurfaceBalance {
                                cell: contact.contact.cell, residual_w, iterations,
                            }
                        }
                    })?;
                (response.surface_temperature_k, response.solid_power_to_bath_w,
                    response.power_to_bath_w, response.conductance_w_k,
                    response.heat_residual_w, response.iterations,
                    response.film.0, response.film.1)
            };
        check()?;
        Ok(FiniteGasContactSample {
            contact: contact.contact,
            specific_enthalpy_j_kg,
            condensed_temperature_k: cell_k,
            condensed_conductivity_w_m_k: conductivity,
            gas_temperature_k: bulk_k,
            film,
            conductance_w_k: conductance,
            heat_into_condensed_w: -power,
            convection: Some(FiniteGasConvectionSample {
                outward_normal: contact.outward_normal,
                surface_temperature_k: surface,
                condensed_power_to_gas_w: solid_power,
                surface_heat_residual_w: residual,
                iterations,
                coefficient,
            }),
        })
    }
}

/// Advance the SAME material/phase/MAC and finite gas with source convection.
///
/// Every epoch retains native k and its phase correction over the half-cell. The
/// common surface root then re-queries gas transport and h at each actual Ts.
/// Positive area and material distance are mandatory; zero gas thickness means
/// convection replaces that film rather than adding a parallel or serial one.
/// Authored unit normals describe stationary geometry, not inferred voxel faces.
///
/// Latent storage, exact solid/mushy occupancy, native momentum, advection and
/// molecular heating retain their original owners. This is a held-cell-phase
/// normal path, not a resolved internal melting front or material k(T) profile.
/// The original explicit timeline applies each signed heat rate and the finite
/// gas's opposite energy once, with combined capacity limits and joint rollback.
/// The returned secant is a heat-rate representation, not a nonlinear Jacobian.
/// No gas momentum, moving geometry, radiation or species transfer is inferred.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_finite_gas_with_material_convection<'t, 'f, F, R, O, T>(
    thermal: &DepositionSimulation<'t>,
    target_s: f64,
    config: StepConfig,
    flow: &MaterialFlowFeedback<'f, F, R, O, T>,
    gas: &IsochoricGasParcel,
    contacts: &[FiniteGasConvectionContact],
    provider: &dyn FiniteGasConvection,
    controls: GasEnergySolve,
    surface: LayerConductionConfig,
) -> Result<
    (DepositionSimulation<'t>, MaterialFlowFeedback<'f, F, R, O, T>,
        IsochoricGasParcel, FiniteGasExchangeReport),
    FiniteGasFlowError<O::Error>,
>
where
    F: Clone,
    T: FlowTransportSnapshot,
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
    O: DepositionObserver,
{
    check_gate(flow.gate, &mut || Ok(()))?;
    if contacts.is_empty() || contacts.len() > flow.input.max_connections
        || surface.max_iterations == 0
        || !surface.absolute_heat_tolerance_w.is_finite()
        || surface.absolute_heat_tolerance_w < 0.0
        || !surface.relative_heat_tolerance.is_finite()
        || surface.relative_heat_tolerance < 0.0
    {
        return Err(invalid("finite gas convection needs bounded contacts and valid surface controls").into());
    }
    let mut scalar = Vec::with_capacity(contacts.len());
    for contact in contacts {
        check_gate(flow.gate, &mut || Ok(()))?;
        let square: f64 = contact.outward_normal.into_iter().map(|x| x * x).sum();
        if !square.is_finite() || (square.sqrt() - 1.0).abs() > 64.0 * f64::EPSILON {
            return Err(invalid("finite gas convection requires authored physical unit normals").into());
        }
        scalar.push(contact.contact);
    }
    advance_finite_gas(thermal, target_s, config, flow, gas, &scalar, controls,
        Some(ConvectionBinding { contacts, provider, config: surface }))
}
