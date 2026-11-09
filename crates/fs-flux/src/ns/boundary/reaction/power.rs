//! Power conjugate to the original admitted data, without boundary requery.

use super::{
    BoundaryReactions, DiscreteBoundaryReaction, StokesBoundaryError, StokesPressureReference,
    Sums, invalid,
};

/// Signed mechanical power ON the fluid, per depth divided by the SAME
/// reference density [m4/s3]. Positive values supply mechanical work.
/// This sign is opposite to the fluid-on-exterior boundary wrenches.
///
/// Prescribed work is the discrete source functional, not a pointwise Cauchy
/// stress estimate: R0*m0(g)+R1*m1(g) plus the full numerical Nitsche flux
/// dotted with the original prescribed trace, or physical Robin effort dotted
/// with the wall trace. The adjoint-consistency terms remain in R. This does
/// not assert that an arbitrary Nitsche discretization conserves energy.
/// Advected kinetic energy, physical slip loss, numerical loss and heat stay
/// separate; neither a moving domain nor whole-scheme passivity is inferred.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundaryPower {
    /// Original prescribed velocity or Navier-slip wall work, including normal
    /// pressure work when the prescribed boundary carries throughflow.
    pub prescribed_velocity: f64,
    /// Literal supplied natural traction dotted with the accepted fluid trace.
    /// Includes any subsequently resolved original pressure-port contribution.
    pub imposed_traction: f64,
    /// Selected dissipative-backflow stress dotted with the actual fluid trace.
    /// Nonpositive, retained separately; it is not automatic material heating.
    pub backflow_resistance: f64,
    /// Selected references: None retains a traction datum; Some adds an
    /// explicit kinematic-pressure offset to the returned mean-zero field.
    pub pressure_offsets: Vec<(usize, Option<f64>)>,
}

impl DiscreteBoundaryReaction {
    /// Actual admitted-boundary work in the solve's velocity frame. No new
    /// motion, source sampling, solve or quadrature is performed. A mean-zero
    /// component requires an explicit offset, including Some(0); a component
    /// with a traction datum requires None. Offsets change prescribed work by
    /// -offset*integral(g.n), retaining the original discrete volume trace.
    ///
    /// # Errors
    /// Refuses incompatible references and unrepresentable retained power.
    pub fn power(
        &self,
        pressure_offset: Option<f64>,
    ) -> Result<BoundaryPower, StokesBoundaryError> {
        let values = self.power_values(pressure_offset)?;
        Ok(BoundaryPower {
            prescribed_velocity: values[0],
            imposed_traction: values[1],
            backflow_resistance: values[2],
            pressure_offsets: vec![(self.component, pressure_offset)],
        })
    }

    fn power_values(&self, pressure_offset: Option<f64>) -> Result<[f64; 3], StokesBoundaryError> {
        let offset = match (self.pressure_reference, pressure_offset) {
            (StokesPressureReference::MeanZero, Some(value)) if value.is_finite() => value,
            (StokesPressureReference::TractionDatum, None) => 0.0,
            _ => {
                return Err(invalid(
                    "explicit mean-zero offset or retained traction datum required",
                ));
            }
        };
        let mut prescribed = Sums::<1>::new();
        prescribed.add([self.power[0]]);
        prescribed.add([-offset * self.power[3]]);
        let values = [prescribed.values()[0], self.power[1], self.power[2]];
        if values.iter().any(|value| !value.is_finite()) || values[2] > 0.0 {
            return Err(invalid("unrepresentable original boundary power"));
        }
        Ok(values)
    }
}

impl BoundaryReactions {
    /// Sum original selected boundary powers in canonical edge order, with
    /// explicit component pressure references and compensated accumulation.
    /// Velocities, constitutive data and actual convective split quadrature
    /// are those retained at the accepted native endpoint. In particular,
    /// stationary prescribed boundaries supply zero wall work even when their
    /// weak fluid trace differs; Navier slip loss remains a separate receipt.
    ///
    /// # Errors
    /// Refuses empty/duplicate/nonboundary selections, wrong component axes,
    /// missing/incompatible pressure datums and nonfinite/unrepresentable power.
    pub fn power(
        &self,
        edges: &[usize],
        pressure_offsets: &[Option<f64>],
    ) -> Result<BoundaryPower, StokesBoundaryError> {
        if pressure_offsets.len() != self.components {
            return Err(invalid(
                "component pressure offsets required for boundary power",
            ));
        }
        let mut sums = Sums::<3>::new();
        let mut referenced = vec![false; self.components];
        for index in self.selected_transport_edges(edges)? {
            let face = &self.edges[index];
            sums.add(face.power_values(pressure_offsets[face.component])?);
            referenced[face.component] = true;
        }
        let values = sums.values();
        if values.iter().any(|value| !value.is_finite()) || values[2] > 0.0 {
            return Err(invalid("unrepresentable boundary power sum"));
        }
        Ok(BoundaryPower {
            prescribed_velocity: values[0],
            imposed_traction: values[1],
            backflow_resistance: values[2],
            pressure_offsets: referenced
                .iter()
                .enumerate()
                .filter_map(|(component, &selected)| {
                    selected.then_some((component, pressure_offsets[component]))
                })
                .collect(),
        })
    }
}
