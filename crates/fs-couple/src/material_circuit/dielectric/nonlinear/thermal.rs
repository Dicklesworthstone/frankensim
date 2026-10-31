//! Retained thermodynamic authority for nonlinear instantaneous storage and
//! affine polarization potentials. The native graph and charge owner are shared.

use super::*;
use fs_material::field_storage::dielectric::{
    reversible::thermal::ThermalReversibleDielectricMaterial, thermal::ThermalDebyeCapacitor,
};

/// Complete original laws for D_inf(E,T)=s(T)D_ref(E) and native Debye poles.
/// This authority is immutable and remains attached to every rebuilt circuit.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermalDebyeStorage {
    polarization: ThermalDebyeCapacitor,
    instantaneous: ThermalReversibleDielectricMaterial,
}

fn close(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 512.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl ThermalDebyeStorage {
    /// Join explicitly admitted original potentials on the SAME temperature
    /// interval. The low-field background and its derivative must agree with
    /// the scaled nonlinear law throughout that interval; it is not a second
    /// instantaneous energy store. Original affine pole laws may vary with T.
    pub fn new(
        polarization: &ThermalDebyeCapacitor,
        instantaneous: &ThermalReversibleDielectricMaterial,
    ) -> Result<Self, MaterialCircuitError> {
        if polarization.card().content_hash() != instantaneous.card().content_hash()
            || polarization.temperature_axis() != instantaneous.temperature_axis()
            || polarization.temperature_range_k() != instantaneous.temperature_range_k()
        {
            return Err(invalid(
                "nonlinear dielectric thermal potentials require the same original card, axis and complete interval",
            ));
        }
        let geometry = polarization.original().area_m2() / polarization.original().gap_m();
        let reference = instantaneous
            .reference()
            .evaluate(0.0)?
            .differential_permittivity_f_m;
        let derivative = finite(finite(reference * geometry)? * instantaneous.scale_slope_per_k())?;
        if !close(polarization.capacitance_slopes_f_k()[0], derivative) {
            return Err(invalid(
                "Debye background derivative contradicts the scaled instantaneous potential",
            ));
        }
        for temperature in polarization.temperature_range_k() {
            let pole_source = polarization.at_temperature(temperature)?;
            let field = instantaneous.at_temperature(temperature)?;
            let actual = finite(
                field
                    .dielectric()
                    .evaluate(0.0)?
                    .differential_permittivity_f_m
                    * geometry,
            )?;
            if !close(actual, pole_source.source().instantaneous_capacitance_f()) {
                return Err(invalid(
                    "Debye background must match the nonlinear tangent over the complete thermal interval",
                ));
            }
        }
        Ok(Self {
            polarization: polarization.clone(),
            instantaneous: instantaneous.clone(),
        })
    }

    pub fn polarization(&self) -> &ThermalDebyeCapacitor {
        &self.polarization
    }
    pub fn instantaneous(&self) -> &ThermalReversibleDielectricMaterial {
        &self.instantaneous
    }

    pub(super) fn admit(
        &self,
        source: &ResolvedDebyeCapacitor,
        instantaneous: &MaterialDielectricCapacitor,
    ) -> Result<(), MaterialCircuitError> {
        super::admit_pair(source, instantaneous)?;
        let poles = self.polarization.at_source(source)?;
        let original = self.instantaneous.at_temperature(poles.temperature_k())?;
        if instantaneous.thermal_material() != Some(&original)
            || instantaneous.material() != original.dielectric()
        {
            return Err(invalid(
                "native instantaneous source must retain its exact original nonlinear thermal potential",
            ));
        }
        Ok(())
    }
}

impl DebyeCircuit {
    /// Assemble the existing native Debye graph with explicitly owned field
    /// entropy and affine polarization potentials. All histories remain native;
    /// the finite thermal owner must receive the same retained authority when
    /// changing temperature. This constructor supplies no external heat sink.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_thermal_reversible_instantaneous(
        graph: CircuitGraph,
        source: ResolvedDebyeCapacitor,
        instantaneous: MaterialDielectricCapacitor,
        thermal_storage: ThermalDebyeStorage,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        polarization_charges_c: Vec<f64>,
        initial_inputs: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        Self::assemble_thermal_instantaneous(
            graph,
            source,
            instantaneous,
            Arc::new(thermal_storage),
            fluxes_wb,
            charges_c,
            polarization_charges_c,
            initial_inputs,
            cancelled,
        )
    }

    // Thermal trials retain the original immutable authority instead of
    // deep-cloning its complete material/pole bundle at every Newton probe.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::material_circuit::dielectric) fn assemble_thermal_instantaneous(
        graph: CircuitGraph,
        source: ResolvedDebyeCapacitor,
        instantaneous: MaterialDielectricCapacitor,
        thermal_storage: Arc<ThermalDebyeStorage>,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        polarization_charges_c: Vec<f64>,
        initial_inputs: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, MaterialCircuitError> {
        poll(&mut cancelled)?;
        thermal_storage.admit(&source, &instantaneous)?;
        let mut circuit = Self::assemble_with_instantaneous(
            graph,
            instantaneous.branch(),
            source,
            Some(instantaneous),
            fluxes_wb,
            charges_c,
            polarization_charges_c,
            initial_inputs,
            &mut cancelled,
        )?;
        Arc::get_mut(&mut circuit.prepared)
            .ok_or_else(|| invalid("unexpected shared nonlinear thermal preparation"))?
            .reversible_thermal = Some(thermal_storage);
        poll(&mut cancelled)?;
        Ok(circuit)
    }

    /// Complete field and polarization thermal authority retained by this graph.
    pub fn thermal_reversible_storage(&self) -> Option<&ThermalDebyeStorage> {
        self.prepared.reversible_thermal.as_deref()
    }
}
