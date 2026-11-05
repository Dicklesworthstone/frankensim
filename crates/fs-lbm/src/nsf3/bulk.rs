//! Source-local molecular bulk viscosity on the original conservative NSF step.
//!
//! For an isotropic Newtonian gas, the trace rate obeys
//! `zeta = (2/3) rho theta0 (1/omega_b - 1/2)` in lattice units.
//! The same admitted rate drives collision and stress work in total-energy
//! transport. No separate heat source, eddy bulk viscosity or molecular
//! relaxation-memory model is introduced.

use super::{NsfCellMaterial, NsfError, NsfGasMedium, NsfGrid3, NsfLaw, THETA0};
use crate::lattice3::Lattice3;
use core::fmt;
use std::sync::Arc;

/// Actual physical coordinates of one native molecular-transport evaluation.
/// Composition and all other material coordinates belong to the retained source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfBulkViscosityPoint {
    /// Local temperature [K].
    pub temperature_k: f64,
    /// Local absolute thermodynamic pressure [Pa], never a lattice gauge.
    pub pressure_pa: f64,
    /// Local physical density [kg/m^3].
    pub density_kg_m3: f64,
    /// Molecular dynamic shear viscosity used by the original flow owner [Pa s].
    /// GasSpec uses its existing transport table; source gases use their exact
    /// local material response. An SGS increment is never included here.
    pub dynamic_viscosity_pa_s: f64,
}

type Evaluate = dyn Fn(NsfBulkViscosityPoint) -> Result<f64, NsfError> + Send + Sync;

/// Retained physical source of instantaneous molecular bulk viscosity [Pa s].
///
/// The callback owns its calibration/material data and must be a deterministic
/// function of the supplied state: predictor queries can repeat and run in
/// parallel. A refusal propagates before any accepted fluid state changes.
/// Clones retain the SAME source through characteristic interval snapshots.
#[derive(Clone)]
pub struct NsfBulkViscositySource {
    citation: Arc<str>,
    evaluate: Arc<Evaluate>,
}

impl fmt::Debug for NsfBulkViscositySource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("NsfBulkViscositySource")
            .field("citation", &self.citation)
            .finish_non_exhaustive()
    }
}

impl NsfBulkViscositySource {
    /// Bind explicit source evidence and its bounded-domain material evaluator.
    /// Positive representable bulk viscosity is admitted at every actual query;
    /// no coefficient is inferred from gas identity or silently clamped.
    pub fn new(
        citation: impl Into<String>,
        evaluate: impl Fn(NsfBulkViscosityPoint) -> Result<f64, NsfError> + Send + Sync + 'static,
    ) -> Result<Self, NsfError> {
        let citation = citation.into();
        if citation.trim().is_empty() {
            return Err(NsfError::Input(
                "bulk viscosity requires an explicit physical source",
            ));
        }
        Ok(Self {
            citation: citation.into(),
            evaluate: Arc::new(evaluate),
        })
    }

    /// Original physical evidence, retained independently of caller buffers.
    pub fn citation(&self) -> &str {
        &self.citation
    }
}

/// Cached material response at the accepted gas state, without another query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfBulkViscositySample {
    /// Actual local molecular state used by the source.
    pub point: NsfBulkViscosityPoint,
    /// Physical coefficient returned by that source [Pa s].
    pub bulk_viscosity_pa_s: f64,
    /// Native trace relaxation rate used by collision and stress work.
    pub trace_relaxation_rate: f64,
}

#[derive(Clone, Copy, Debug)]
pub(super) struct BulkTransport {
    viscosity_pa_s: f64,
    pub(super) trace_rate: f64,
}

impl NsfLaw {
    fn bulk_point(&self, rho: f64, theta: f64, mu_lattice: f64) -> NsfBulkViscosityPoint {
        NsfBulkViscosityPoint {
            temperature_k: self.kelvin_of_theta(theta),
            pressure_pa: self.pressure_pa(rho, theta),
            density_kg_m3: rho * self.rho_ref_kg_m3,
            dynamic_viscosity_pa_s: mu_lattice / self.transport_scale(),
        }
    }

    fn sample_bulk(&self, rho: f64, theta: f64, mu: f64) -> Result<BulkTransport, NsfError> {
        let point = self.bulk_point(rho, theta, mu);
        if ![
            point.temperature_k,
            point.pressure_pa,
            point.density_kg_m3,
            point.dynamic_viscosity_pa_s,
        ]
        .into_iter()
        .all(|v| v.is_finite() && v > 0.0)
        {
            return Err(NsfError::Input(
                "unrepresentable physical bulk-viscosity query",
            ));
        }
        let source = self
            .bulk_source
            .as_ref()
            .expect("explicit bulk source selected");
        let viscosity_pa_s = (source.evaluate)(point)?;
        let viscosity_lattice = viscosity_pa_s * self.transport_scale();
        let trace_rate = 1.0 / (0.5 + 1.5 * viscosity_lattice / (rho * THETA0));
        if !(viscosity_pa_s.is_finite()
            && viscosity_pa_s > 0.0
            && viscosity_lattice.is_finite()
            && viscosity_lattice > 0.0
            && trace_rate.is_finite()
            && trace_rate > 0.0
            && trace_rate < 2.0)
        {
            return Err(NsfError::Input(
                "physical bulk viscosity or trace relaxation is unrepresentable",
            ));
        }
        Ok(BulkTransport {
            viscosity_pa_s,
            trace_rate,
        })
    }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Attach a source-local `zeta(T,p,rho,mu)` before the first native step.
    /// The original EOS, shear transport, thermal conduction, boundaries and
    /// total-energy integrator stay authoritative. Current flux history is
    /// rebuilt with the same trace rate that will drive collision. The source
    /// is retained through original characteristic trial/accept transactions.
    ///
    /// Positive representability and the existing molecular/Mach/temperature
    /// bounds remain required. Their admission is not a stability qualification
    /// of every new bulk law. This instantaneous constitutive model supplies
    /// no frequency-dependent relaxation, extra heat, or SGS bulk stress.
    pub fn with_bulk_viscosity(mut self, source: NsfBulkViscositySource) -> Result<Self, NsfError> {
        if self.steps != 0 || self.law.bulk_source.is_some() {
            return Err(NsfError::Input(
                "attach the bulk-viscosity source once, before stepping",
            ));
        }
        self.law.bulk_source = Some(source);
        let mu = self
            .law
            .source_reference_transport
            .map_or_else(|| self.law.mu_lattice(THETA0), |t| t[0]);
        self.law.source_reference_bulk_pa_s =
            Some(self.law.sample_bulk(1.0, THETA0, mu)?.viscosity_pa_s);
        let (rho, momentum) = self.macro_fields(&self.f);
        let velocity = self.velocity_field(&rho, &momentum, &self.force);
        let material = self.material_field(&rho, &self.theta)?;
        let stress = self.stress_field(
            &self.f,
            &self.force,
            &rho,
            &self.theta,
            &material,
            &self.sgs_transport,
        );
        let (div, wall, open) = self.flux_divergence(
            &self.e,
            &rho,
            &velocity,
            &self.theta,
            &stress,
            &self.sgs_transport,
            self.open_faces.as_ref(),
        )?;
        self.material = material;
        self.div = div;
        self.wall_out = wall;
        self.open_out = open;
        Ok(self)
    }

    /// Source evidence of the explicitly selected local bulk law, if present.
    pub fn bulk_viscosity_source(&self) -> Option<&str> {
        self.law
            .bulk_source
            .as_ref()
            .map(NsfBulkViscositySource::citation)
    }

    /// Accepted local bulk response. `None` denotes the unchanged ratio path.
    pub fn bulk_viscosity_at(
        &self,
        x: [usize; 3],
    ) -> Result<Option<NsfBulkViscositySample>, NsfError> {
        if x.iter().zip(self.spec.dims).any(|(&i, n)| i >= n) {
            return Err(NsfError::Input("bulk-viscosity cell outside the grid"));
        }
        let i = self.index(x);
        if self.solid[i] {
            return Err(NsfError::Input("solid cell has no gas bulk viscosity"));
        }
        let Some(material) = self.material.get(i).and_then(Option::as_ref) else {
            return Ok(None);
        };
        let Some(bulk) = material.bulk else {
            return Ok(None);
        };
        let (rho, _) = super::raw_moments::<L>(&self.f[i]);
        Ok(Some(NsfBulkViscositySample {
            point: self.law.bulk_point(rho, self.theta[i], material.mu_lattice),
            bulk_viscosity_pa_s: bulk.viscosity_pa_s,
            trace_relaxation_rate: bulk.trace_rate,
        }))
    }

    pub(super) fn material_with_bulk(
        &self,
        rho: f64,
        theta: f64,
        material: Option<NsfCellMaterial>,
    ) -> Result<Option<NsfCellMaterial>, NsfError> {
        if self.law.bulk_source.is_none() {
            return Ok(material);
        }
        let mut material = if matches!(&self.spec.gas, NsfGasMedium::CaloricallyPerfect(_)) {
            // Keep the original GasSpec table and exact cv*theta arithmetic.
            // The original gas evaluator supplies the thermodynamic snapshot;
            // it is not substituted for the native transport table.
            self.spec.gas.resolve(&self.law, rho, theta)?;
            NsfCellMaterial {
                state: self.spec.gas.at(
                    self.law.kelvin_of_theta(theta),
                    self.law.pressure_pa(rho, theta),
                )?,
                mu_lattice: self.law.mu_lattice(theta),
                sensible_energy_lattice: (1.0 / (self.law.gamma - 1.0)) * theta,
                bulk: None,
            }
        } else {
            material.ok_or(NsfError::Input(
                "source gas is missing its native material snapshot",
            ))?
        };
        material.bulk = Some(self.law.sample_bulk(rho, theta, material.mu_lattice)?);
        Ok(Some(material))
    }
}
