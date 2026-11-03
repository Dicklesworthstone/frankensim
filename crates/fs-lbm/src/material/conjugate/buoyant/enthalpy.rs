//! Original caloric laws on the existing material/flow fixed point and clock.

mod transient;
#[cfg(feature = "enthalpy-coupling")]
pub use transient::MaterialViscousHeatMarch;
pub use transient::{
    MaterialEnthalpyMarch, MaterialEnthalpySolidMarch, MaterialEnthalpyStepReport,
    MaterialEnthalpyTransientConfig,
};

use super::{
    MaterialBuoyancyConfig, MaterialBuoyancyError, MaterialBuoyancySolution, MaterialChtFluid,
    SolidConductivityLaw, poll,
};
use crate::conjugate::{
    AdvectiveEnthalpyField, AffineSpecificEnthalpy, ConvectionScheme, EnergyConfig, FlowField, SimpleConfig, SolidMaterial,
    ThermalFace, ThermalSetup, VoxelDomain,
};
use crate::material::transport::boussinesq::BoussinesqConfig;
use fs_exec::CancelGate;
use fs_material::{
    MaterialError,
    fluid_transport::{FluidTransport, FluidTransportState},
    gas_species::input::temperature::IsobaricGasLaw,
    phase::ResolvedHeatCapacityCurve,
    thermochemical_gas::energy::enthalpy::GasEnthalpyLaw,
};
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityLaw;
use fs_qty::{Pressure, Temperature};

/// One immutable, single-phase caloric owner at fixed pressure/composition.
///
/// The returned h is `h(temperature_k) - h(reference_temperature_k)`, evaluated
/// directly from the original law, and Cp is that same law's positive tangent.
/// An arbitrary formation-energy datum must not enter the numerical flux.
/// The transport resolver must describe this SAME fluid and source context;
/// this interface cannot infer a material identity from numerical coefficients.
/// That source/context correspondence is an explicit caller contract. The
/// heat-capacity chart's chord derivative can differ from the raw source Cp
/// within its declared interpolation approximation; no arbitrary pointwise
/// equality between those two capacities is imposed.
pub trait ConvectionCaloricLaw {
    /// Evaluate original-law relative enthalpy and its local positive Cp.
    /// Both temperatures must remain inside the same admitted source branch.
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError>;
}

fn source_error(error: impl core::fmt::Display) -> MaterialError {
    MaterialError::State {
        what: format!("convection caloric law: {error}"),
    }
}

impl ConvectionCaloricLaw for ResolvedHeatCapacityCurve {
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError> {
        Ok(AffineSpecificEnthalpy {
            temperature_k,
            enthalpy_j_kg: self
                .specific_enthalpy_change(reference_temperature_k, temperature_k)
                .map_err(source_error)?,
            heat_capacity_j_kg_k: self
                .secant_heat_capacity(temperature_k, temperature_k)
                .map_err(source_error)?,
        })
    }
}

impl ConvectionCaloricLaw for IsobaricHeatCapacityLaw {
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError> {
        Ok(AffineSpecificEnthalpy {
            temperature_k,
            enthalpy_j_kg: self
                .specific_enthalpy_change_j_kg(reference_temperature_k, temperature_k)
                .map_err(source_error)?,
            heat_capacity_j_kg_k: self
                .specific_heat_cp_j_kg_k(temperature_k)
                .map_err(source_error)?,
        })
    }
}

impl ConvectionCaloricLaw for IsobaricGasLaw {
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError> {
        let temperature = Temperature::new(temperature_k);
        let endpoint = self.at_temperature(temperature)?;
        Ok(AffineSpecificEnthalpy {
            temperature_k,
            enthalpy_j_kg: self.specific_enthalpy_change(
                Temperature::new(reference_temperature_k), temperature,
            )?,
            heat_capacity_j_kg_k: endpoint.state().specific_heat_cp,
        })
    }
}

impl ConvectionCaloricLaw for GasEnthalpyLaw {
    fn tangent_at_temperature(
        &self,
        reference_temperature_k: f64,
        temperature_k: f64,
    ) -> Result<AffineSpecificEnthalpy, MaterialError> {
        // These calls retain this owner's explicit caloric bracket. Their
        // absolute h values are never subtracted to obtain the transported h.
        self.enthalpy_at_temperature(reference_temperature_k)
            .map_err(source_error)?;
        self.enthalpy_at_temperature(temperature_k)
            .map_err(source_error)?;
        let pressure = self
            .context()
            .axes()
            .get("pressure")
            .copied()
            .ok_or_else(|| source_error("gas caloric context has no pressure"))?;
        let properties = self
            .model()
            .thermodynamics()
            .evaluate_properties(Temperature::new(temperature_k), Pressure::new(pressure))
            .map_err(source_error)?;
        let change = self
            .model()
            .thermodynamics()
            .specific_enthalpy_change(
                Temperature::new(reference_temperature_k),
                Temperature::new(temperature_k),
            )
            .map_err(source_error)?;
        Ok(AffineSpecificEnthalpy {
            temperature_k,
            enthalpy_j_kg: change.value(),
            heat_capacity_j_kg_k: properties.mass_specific_cp().value(),
        })
    }
}

/// Caloric coordinates and nonlinear acceptance, independent of Krylov tolerance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialEnthalpyConfig {
    /// Supported common reference for `h(T)-h(T_reference)`, K. This is a
    /// numerical datum, not the buoyancy reference or an ambient reservoir.
    pub reference_temperature_k: f64,
    /// Absolute maximum defect [W] of EACH advective face and EACH cell's
    /// net advective heat after re-evaluating the true endpoint caloric law.
    pub maximum_flux_defect_w: f64,
}

/// Defects of the actual caloric law against the solved affine face operator.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MaterialEnthalpyReport {
    /// Largest absolute true-minus-affine advective face heat, W.
    pub maximum_face_flux_defect_w: f64,
    /// Largest absolute net true-minus-affine advective cell heat, W.
    pub maximum_cell_flux_defect_w: f64,
}

/// Applied and re-evaluated caloric fields on the same material fixed point.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialEnthalpyState {
    /// None only for solid voxels; those keep the existing conduction owner.
    pub applied: Vec<Option<AffineSpecificEnthalpy>>,
    /// True caloric tangents re-evaluated at the returned temperatures.
    pub endpoint: Vec<Option<AffineSpecificEnthalpy>>,
    /// Source-law h at every actual inlet/backflow temperature, on the same datum.
    pub boundary_j_kg: [Option<f64>; 6],
    /// Local nonlinear heat defects admitted at the returned fixed point.
    pub report: MaterialEnthalpyReport,
}

/// Original flow/solid/fluid result plus the borrowed original caloric owner.
pub struct MaterialEnthalpySolution<'a, F, C: ConvectionCaloricLaw> {
    /// Original solved flow, heat, fluid sources and optional solid sources.
    pub convection: MaterialBuoyancySolution<F>,
    /// Same immutable caloric owner used by every cell and boundary.
    pub caloric: &'a C,
    /// Common caloric reference and local heat-flux acceptance tolerance.
    pub config: MaterialEnthalpyConfig,
    /// The energy report audits the applied affine law. These local defects
    /// separately bound its mismatch to the true endpoint caloric law; their
    /// cancellation in a global balance is never used to accept a solution.
    pub enthalpy: MaterialEnthalpyState,
}

#[derive(Clone, Copy)]
pub(super) struct CaloricCoupling<'a> {
    pub law: &'a dyn ConvectionCaloricLaw,
    pub config: MaterialEnthalpyConfig,
}

impl CaloricCoupling<'_> {
    pub(super) fn validate(self, energy: &EnergyConfig) -> Result<(), MaterialBuoyancyError> {
        if !self.config.reference_temperature_k.is_finite()
            || self.config.reference_temperature_k <= 0.0
            || !self.config.maximum_flux_defect_w.is_finite()
            || self.config.maximum_flux_defect_w < 0.0
            || !matches!(energy.scheme, ConvectionScheme::Upwind | ConvectionScheme::EnthalpyMinmod)
        {
            return Err(MaterialBuoyancyError::Invalid(
                "source enthalpy requires positive datum temperature, finite nonnegative flux tolerance and upwind or limited enthalpy transport",
            ));
        }
        Ok(())
    }

    fn sample(
        self,
        temperature_k: f64,
        reference: &FluidTransportState,
        controls: BoussinesqConfig<3>,
    ) -> Result<AffineSpecificEnthalpy, MaterialBuoyancyError> {
        let sample = self
            .law
            .tangent_at_temperature(self.config.reference_temperature_k, temperature_k)
            .map_err(|source| MaterialBuoyancyError::Caloric {
                temperature_k,
                source,
            })?;
        let cp = sample.heat_capacity_j_kg_k;
        let reference_cp = reference.specific_heat_j_kg_k();
        if sample.temperature_k != temperature_k
            || !temperature_k.is_finite()
            || temperature_k <= 0.0
            || !sample.enthalpy_j_kg.is_finite()
            || !cp.is_finite()
            || cp <= 0.0
            || ((cp - reference_cp) / reference_cp).abs()
                > controls.max_relative_specific_heat_departure
        {
            return Err(MaterialBuoyancyError::Invalid(
                "caloric tangent must match its temperature and original positive-Cp departure bound",
            ));
        }
        Ok(sample)
    }

    pub(super) fn field(
        self,
        domain: &VoxelDomain,
        temperatures: &[f64],
        reference: &FluidTransportState,
        controls: BoussinesqConfig<3>,
        gate: &CancelGate,
    ) -> Result<Vec<Option<AffineSpecificEnthalpy>>, MaterialBuoyancyError> {
        let mut field = Vec::with_capacity(domain.cell_count());
        self.write_field(domain, temperatures, reference, controls, &mut field, gate)?;
        Ok(field)
    }

    fn write_field(
        self,
        domain: &VoxelDomain,
        temperatures: &[f64],
        reference: &FluidTransportState,
        controls: BoussinesqConfig<3>,
        field: &mut Vec<Option<AffineSpecificEnthalpy>>,
        gate: &CancelGate,
    ) -> Result<(), MaterialBuoyancyError> {
        field.clear();
        for (cell, &temperature) in temperatures.iter().enumerate() {
            poll(gate)?;
            field.push(if domain.is_fluid(cell) {
                Some(self.sample(temperature, reference, controls)?)
            } else {
                None
            });
        }
        Ok(())
    }

    pub(super) fn boundaries(
        self,
        setup: &ThermalSetup,
        reference: &FluidTransportState,
        controls: BoussinesqConfig<3>,
        gate: &CancelGate,
    ) -> Result<[Option<f64>; 6], MaterialBuoyancyError> {
        let mut values = [None; 6];
        for (face, rule) in setup.faces.iter().enumerate() {
            poll(gate)?;
            if let ThermalFace::Inflow { temperature }
            | ThermalFace::Outflow {
                backflow_temperature: temperature,
            } = *rule
            {
                values[face] = Some(self.sample(temperature, reference, controls)?.enthalpy_j_kg);
            }
        }
        Ok(values)
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn defects(
    domain: &VoxelDomain,
    flow: &FlowField,
    density: f64,
    scheme: ConvectionScheme,
    temperatures: &[f64],
    applied: &[Option<AffineSpecificEnthalpy>],
    endpoint: &[Option<AffineSpecificEnthalpy>],
    gate: &CancelGate,
) -> Result<MaterialEnthalpyReport, MaterialBuoyancyError> {
    let mut cell_defects = vec![0.0; domain.cell_count()];
    defect_field(domain, flow, density, scheme, temperatures, applied, endpoint, &mut cell_defects, gate)
}

#[allow(clippy::too_many_arguments)]
fn defect_field(
    domain: &VoxelDomain,
    flow: &FlowField,
    density: f64,
    scheme: ConvectionScheme,
    temperatures: &[f64],
    applied: &[Option<AffineSpecificEnthalpy>],
    endpoint: &[Option<AffineSpecificEnthalpy>],
    cell_defects: &mut [f64],
    gate: &CancelGate,
) -> Result<MaterialEnthalpyReport, MaterialBuoyancyError> {
    let mut report = MaterialEnthalpyReport::default();
    cell_defects.fill(0.0);
    // The original source owner re-evaluates BOTH the true caloric law and
    // its limited traces. A converged tangent alone cannot accept a stale
    // nonlocal reconstruction. Each image is reconstructed once per trial.
    let offsets = if scheme == ConvectionScheme::EnthalpyMinmod {
        Some((
            AdvectiveEnthalpyField { cells: applied, boundary_j_kg: [None; 6] }
                .limited_offsets(domain, gate)?,
            AdvectiveEnthalpyField { cells: endpoint, boundary_j_kg: [None; 6] }
                .limited_offsets(domain, gate)?,
        ))
    } else {
        None
    };
    for cell in 0..domain.cell_count() {
        poll(gate)?;
        for face in 0..6 {
            let neighbor = domain.neighbor(cell, face);
            if neighbor.is_some_and(|next| next < cell) {
                continue;
            }
            let volume_rate = flow.outward(domain, cell, face);
            if volume_rate == 0.0 {
                continue;
            }
            let donor = if volume_rate > 0.0 {
                cell
            } else if let Some(next) = neighbor {
                next
            } else {
                continue;
            }; // Prescribed incoming h is already exact.
            let previous = applied[donor].expect("admitted fluid donor");
            let actual = endpoint[donor].expect("admitted fluid endpoint");
            // Subtract the common enthalpy datum before applying the sensible
            // increment. Adding that increment to a large absolute h first
            // can round it away and falsely certify a zero caloric defect.
            let mut specific_difference = (-previous.heat_capacity_j_kg_k).mul_add(
                temperatures[donor] - previous.temperature_k,
                actual.enthalpy_j_kg - previous.enthalpy_j_kg,
            );
            if neighbor.is_some() && let Some((applied_offsets, endpoint_offsets)) = &offsets {
                let donor_face = if volume_rate > 0.0 { face } else { face ^ 1 };
                let sign = if donor_face % 2 == 0 { -1.0 } else { 1.0 };
                specific_difference += sign
                    * (endpoint_offsets[donor][face / 2] - applied_offsets[donor][face / 2]);
            }
            let difference = density * volume_rate * specific_difference;
            cell_defects[cell] += difference;
            if let Some(next) = neighbor {
                cell_defects[next] -= difference;
            }
            if !difference.is_finite()
                || !cell_defects[cell].is_finite()
                || neighbor.is_some_and(|next| !cell_defects[next].is_finite())
            {
                return Err(MaterialBuoyancyError::Invalid(
                    "nonlinear enthalpy face defect is unrepresentable",
                ));
            }
            report.maximum_face_flux_defect_w =
                report.maximum_face_flux_defect_w.max(difference.abs());
        }
    }
    report.maximum_cell_flux_defect_w = cell_defects.iter().copied().map(f64::abs).fold(0.0, f64::max);
    Ok(report)
}

impl<F: FluidTransport + ?Sized> MaterialChtFluid<'_, F> {
    /// Steady natural/mixed convection with the original fluid's full h(T).
    ///
    /// The existing material fixed point still owns rho(T), viscosity, k and
    /// buoyancy; constant reference density still owns continuity/inertia.
    /// Both original Boussinesq density AND Cp departure bounds remain enforced,
    /// including the caloric tangent. This opt-in entry does not extend those
    /// physical approximation bounds or infer compressibility/pressure work.
    ///
    /// `caloric` and `resolve` must bind the same fluid and fixed source context.
    /// Each open face queries that same law at its declared inlet/backflow
    /// temperature. The borrowed source owner and both caloric epochs survive
    /// in the returned result. Positive single-phase Cp and conservative
    /// source-enthalpy transport are required; phase plateaux retain their
    /// separate h-state owner.
    /// The original transport resolver also admits every incident fluid cell
    /// at those boundary temperatures, including source support and both the
    /// reference-density and source-Cp departure bounds before solving flow.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_natural_convection_with_enthalpy<
        'a,
        S: FluidTransport,
        C: ConvectionCaloricLaw,
    >(
        &self,
        domain: &VoxelDomain,
        solids: &[SolidMaterial],
        flow: &SimpleConfig,
        thermal: &ThermalSetup,
        energy: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: MaterialBuoyancyConfig,
        caloric: &'a C,
        enthalpy: MaterialEnthalpyConfig,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialEnthalpySolution<'a, S, C>, MaterialBuoyancyError> {
        let (convection, state) = self.solve_natural_convection_inner(
            domain,
            solids,
            None,
            flow,
            thermal,
            energy,
            initial_temperature_k,
            config,
            Some(CaloricCoupling {
                law: caloric,
                config: enthalpy,
            }),
            gate,
            resolve,
        )?;
        Ok(MaterialEnthalpySolution {
            convection,
            caloric,
            config: enthalpy,
            enthalpy: state.expect("requested caloric coupling"),
        })
    }

    /// The same source-caloric flow with original solid conductivity sources.
    /// Solids retain the existing axis-aligned half-cell/contact operator and
    /// endpoint source checks; they never receive advective fluid enthalpy.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_natural_convection_with_material_solids_and_enthalpy<
        'a,
        S: FluidTransport,
        C: ConvectionCaloricLaw,
    >(
        &self,
        domain: &VoxelDomain,
        sources: &[SolidConductivityLaw<'_>],
        flow: &SimpleConfig,
        thermal: &ThermalSetup,
        energy: &EnergyConfig,
        initial_temperature_k: &[f64],
        config: MaterialBuoyancyConfig,
        caloric: &'a C,
        enthalpy: MaterialEnthalpyConfig,
        gate: &CancelGate,
        resolve: impl FnMut(usize, f64, &FluidTransportState) -> Result<S, MaterialError>,
    ) -> Result<MaterialEnthalpySolution<'a, S, C>, MaterialBuoyancyError> {
        poll(gate)?;
        let solids: Vec<_> = sources
            .iter()
            .map(|source| {
                let k = source.reference_conductivity_w_m_k();
                SolidMaterial::new("source-bound solid", k[0]).with_orthotropic(k)
            })
            .collect();
        let (convection, state) = self.solve_natural_convection_inner(
            domain,
            &solids,
            Some(sources),
            flow,
            thermal,
            energy,
            initial_temperature_k,
            config,
            Some(CaloricCoupling {
                law: caloric,
                config: enthalpy,
            }),
            gate,
            resolve,
        )?;
        Ok(MaterialEnthalpySolution {
            convection,
            caloric,
            config: enthalpy,
            enthalpy: state.expect("requested caloric coupling"),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::FlowFace;

    #[test]
    fn source_enthalpy_interior_face_defect_cannot_hide_in_zero_global_boundary_defect() {
        let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[0] = FlowFace::Fixed;
        faces[1] = FlowFace::Fixed;
        let flow = FlowField::from_face_velocity(
            &domain,
            faces,
            |axis, _| if axis == 0 { 0.1 } else { 0.0 },
        )
        .unwrap();
        let tangent = AffineSpecificEnthalpy {
            temperature_k: 300.0,
            enthalpy_j_kg: 0.0,
            heat_capacity_j_kg_k: 2.0,
        };
        let applied = [Some(tangent); 3];
        let mut endpoint = applied;
        // Only an interior donor differs. Both external boundaries carry their
        // exact original h, and summing the two affected rows gives zero.
        endpoint[1].as_mut().unwrap().enthalpy_j_kg = 4.0;
        let report = defects(
            &domain,
            &flow,
            1.0,
            ConvectionScheme::Upwind,
            &[300.0; 3],
            &applied,
            &endpoint,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(report.maximum_face_flux_defect_w, 0.4);
        assert_eq!(report.maximum_cell_flux_defect_w, 0.4);
        assert_eq!(applied[0], endpoint[0]);
        assert_eq!(applied[2], endpoint[2]);
    }

    #[test]
    fn limited_trace_defect_is_checked_even_when_every_caloric_tangent_is_exact() {
        let domain = VoxelDomain::new(5, 1, 1, 1.0).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[0] = FlowFace::Fixed;
        faces[1] = FlowFace::Fixed;
        let flow = FlowField::from_face_velocity(
            &domain, faces, |axis, _| if axis == 0 { 0.1 } else { 0.0 },
        ).unwrap();
        let tangent = |temperature_k| Some(AffineSpecificEnthalpy {
            temperature_k,
            enthalpy_j_kg: 2.0 * (temperature_k - 300.0),
            heat_capacity_j_kg_k: 2.0,
        });
        let applied = [300.0; 5].map(tangent);
        let temperatures = [300.0, 301.0, 302.0, 303.0, 304.0];
        let endpoint = temperatures.map(tangent);
        let gate = CancelGate::new();
        // The exact affine caloric law has zero donor-value mismatch. Its old
        // uniform reconstruction is nevertheless missing 1 J/kg on each
        // interior outgoing face of this resolved linear enthalpy profile.
        let donor = defects(&domain, &flow, 1.0, ConvectionScheme::Upwind,
            &temperatures, &applied, &endpoint, &gate).unwrap();
        assert_eq!(donor, MaterialEnthalpyReport::default());
        let limited = defects(&domain, &flow, 1.0, ConvectionScheme::EnthalpyMinmod,
            &temperatures, &applied, &endpoint, &gate).unwrap();
        assert!((limited.maximum_face_flux_defect_w - 0.1).abs() < 1e-15);
        assert!((limited.maximum_cell_flux_defect_w - 0.1).abs() < 1e-15);
        let converged = defects(&domain, &flow, 1.0, ConvectionScheme::EnthalpyMinmod,
            &temperatures, &endpoint, &endpoint, &gate).unwrap();
        assert_eq!(converged, MaterialEnthalpyReport::default());
    }

    #[test]
    fn small_caloric_defects_survive_large_positive_and_negative_enthalpy_datums() {
        let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[0] = FlowFace::Fixed;
        faces[1] = FlowFace::Fixed;
        let temperature_increment = 2_f64.powi(-16);
        let temperatures = [300.0, 300.0 + temperature_increment, 300.0];
        for reverse in [false, true] {
            let flow = FlowField::from_face_velocity(&domain, faces, |axis, _| {
                if axis == 0 {
                    if reverse { -0.125 } else { 0.125 }
                } else {
                    0.0
                }
            })
            .unwrap();
            for datum in [2_f64.powi(40), -2_f64.powi(40)] {
                let tangent = AffineSpecificEnthalpy {
                    temperature_k: 300.0,
                    enthalpy_j_kg: datum,
                    heat_capacity_j_kg_k: 2.0,
                };
                let applied = [Some(tangent); 3];
                let mut endpoint = applied;
                endpoint[1].as_mut().unwrap().temperature_k = temperatures[1];
                // The represented endpoint h is unchanged. A tangent with a
                // nonzero Cp nevertheless predicts a small positive increment.
                // Evaluating that tangent in absolute h would hide the defect.
                assert_eq!(tangent.at_temperature(temperatures[1]), datum);
                let expected = 2.0 * 0.125 * 2.0 * temperature_increment;
                for scheme in [ConvectionScheme::Upwind, ConvectionScheme::EnthalpyMinmod] {
                    let report = defects(
                        &domain,
                        &flow,
                        2.0,
                        scheme,
                        &temperatures,
                        &applied,
                        &endpoint,
                        &CancelGate::new(),
                    )
                    .unwrap();
                    assert_eq!(report.maximum_face_flux_defect_w, expected);
                    assert_eq!(report.maximum_cell_flux_defect_w, expected);
                }
            }
        }
    }
}
