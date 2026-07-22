//! Thermal boundary rows and the boundary partition they live on.
//!
//! # The three families
//!
//! | family | prescribes | units | weak-form contribution |
//! | --- | --- | --- | --- |
//! | Dirichlet | `T` | K | eliminated (row/column lift) |
//! | Neumann | `q_n = (−k∇T)·n` | W/m² | `− ∫_Γ v q_n dA` on the load |
//! | Robin | `h`, `T_ref` | W/(m²·K), K | `+ ∫_Γ h v T dA` on the operator, `+ ∫_Γ h v T_ref dA` on the load |
//!
//! **Sign convention, declared once:** `q_n` is the OUTWARD heat flux
//! density — positive means heat LEAVING the domain. `fs-scenario`'s
//! `Physics::Thermal / BcKind::Neumann` row fixes the dimensions
//! (W/m²) but not the sign, so this is the crate's own declaration and
//! the lowering in [`ThermalBc::from_scenario_row`] applies it verbatim.
//!
//! # The Robin seam
//!
//! `fs-scenario`'s Robin row carries the transfer coefficient ONLY:
//! `expectation(Physics::Thermal, BcKind::Robin) == Value(HTC)`. There
//! is no companion `T_ref` field on `BoundaryCondition`. This crate
//! therefore requires `T_ref` to be named at the lowering call — either
//! explicitly, or from the scenario's
//! `Environment::ambient_temperature`. That is the seam the E05
//! correlation rung plugs into: a correlation computes `h`, the
//! reference temperature stays a declared property of the row, and the
//! solve never has to guess which ambient a coefficient was fitted
//! against.

use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Expectation, Physics, expectation};
use fs_scenario::scenario::Environment;

use crate::field::ScalarField;
use crate::mesh::{BoundaryFace, ConductionMesh};
use crate::{ConductionError, HEAT_FLUX_DIMS, HTC_DIMS, TEMPERATURE_DIMS};

/// One thermal boundary condition.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalBc {
    /// Prescribed temperature, K.
    Dirichlet {
        /// Temperature field, K.
        temperature: ScalarField,
    },
    /// Prescribed OUTWARD heat flux density, W/m² (positive = leaving).
    Neumann {
        /// Outward flux field, W/m².
        outward_flux: ScalarField,
    },
    /// Convective transfer: `(−k∇T)·n = h (T − T_ref)`.
    Robin {
        /// Transfer coefficient, W/(m²·K); must be positive.
        htc: ScalarField,
        /// Reference (ambient) temperature, K.
        t_ref: ScalarField,
    },
}

impl ThermalBc {
    /// A stable tag for structured logs.
    #[must_use]
    pub const fn tag(&self) -> &'static str {
        match self {
            ThermalBc::Dirichlet { .. } => "dirichlet",
            ThermalBc::Neumann { .. } => "neumann",
            ThermalBc::Robin { .. } => "robin",
        }
    }

    /// A uniform Dirichlet temperature.
    ///
    /// # Errors
    /// [`ConductionError::NonFinite`].
    pub fn dirichlet(temperature_k: f64) -> Result<ThermalBc, ConductionError> {
        Ok(ThermalBc::Dirichlet {
            temperature: ScalarField::uniform("dirichlet temperature", temperature_k)?,
        })
    }

    /// A uniform outward heat flux.
    ///
    /// # Errors
    /// [`ConductionError::NonFinite`].
    pub fn neumann(outward_flux_w_m2: f64) -> Result<ThermalBc, ConductionError> {
        Ok(ThermalBc::Neumann {
            outward_flux: ScalarField::uniform("neumann outward flux", outward_flux_w_m2)?,
        })
    }

    /// An adiabatic (zero-flux) row. The same thing the natural boundary
    /// condition gives, spelled out so a model never relies on silence.
    #[must_use]
    pub fn adiabatic() -> ThermalBc {
        ThermalBc::Neumann {
            outward_flux: ScalarField::Uniform(0.0),
        }
    }

    /// A uniform convective row.
    ///
    /// # Errors
    /// [`ConductionError::NonFinite`] or
    /// [`ConductionError::Config`] for a non-positive `h`.
    pub fn robin(htc_w_m2k: f64, t_ref_k: f64) -> Result<ThermalBc, ConductionError> {
        if !(htc_w_m2k.is_finite() && htc_w_m2k > 0.0) {
            return Err(ConductionError::Config {
                parameter: "robin htc",
                what: format!("h = {htc_w_m2k} W/(m^2 K) must be finite and positive"),
            });
        }
        Ok(ThermalBc::Robin {
            htc: ScalarField::uniform("robin htc", htc_w_m2k)?,
            t_ref: ScalarField::uniform("robin reference temperature", t_ref_k)?,
        })
    }

    /// Validate against a mesh's vertex count and the physical domain of
    /// each quantity.
    ///
    /// # Errors
    /// [`ConductionError::FieldLength`], [`ConductionError::NonFinite`],
    /// or [`ConductionError::Config`] for a non-positive transfer
    /// coefficient.
    pub fn validate(&self, vertex_count: usize) -> Result<(), ConductionError> {
        match self {
            ThermalBc::Dirichlet { temperature } => {
                temperature.validate("dirichlet temperature", vertex_count)
            }
            ThermalBc::Neumann { outward_flux } => {
                outward_flux.validate("neumann outward flux", vertex_count)
            }
            ThermalBc::Robin { htc, t_ref } => {
                htc.validate("robin htc", vertex_count)?;
                t_ref.validate("robin reference temperature", vertex_count)?;
                let bad = match htc {
                    ScalarField::Uniform(v) => (*v <= 0.0).then_some(*v),
                    ScalarField::Nodal(values) => values.iter().copied().find(|v| *v <= 0.0),
                };
                if let Some(v) = bad {
                    return Err(ConductionError::Config {
                        parameter: "robin htc",
                        what: format!("h = {v} W/(m^2 K) must be positive everywhere"),
                    });
                }
                Ok(())
            }
        }
    }

    /// Lower an `fs-scenario` thermal row into a solvable condition.
    ///
    /// `t_ref` is REQUIRED for Robin rows and ignored otherwise: the
    /// scenario row carries only `h` (see the module docs). Passing
    /// `None` falls back to `environment.ambient_temperature`, which is
    /// a declaration, not a default — the caller chose to say "this
    /// coefficient is referenced to the scenario ambient".
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a non-thermal physics, a
    /// kind outside this crate's vocabulary, a missing value, or a
    /// time-signal/profile/typed payload (a steady solve consumes no
    /// time histories); [`ConductionError::Dimensions`] when the value's
    /// SI exponents disagree with `fs-scenario`'s own contract.
    pub fn from_scenario_row(
        row: &BoundaryCondition,
        environment: &Environment,
        t_ref: Option<f64>,
    ) -> Result<ThermalBc, ConductionError> {
        if row.physics != Physics::Thermal {
            return Err(ConductionError::ScenarioRow {
                region: row.region.clone(),
                what: format!("physics is {:?}, not Thermal", row.physics),
                fix: "attach the row to Physics::Thermal".to_string(),
            });
        }
        let expected = match expectation(row.physics, row.kind) {
            Expectation::Value(dims) => dims,
            other => {
                return Err(ConductionError::ScenarioRow {
                    region: row.region.clone(),
                    what: format!(
                        "fs-scenario expectation for {:?} is {other:?}, which this crate \
                         does not solve",
                        row.kind
                    ),
                    fix: "use Dirichlet, Neumann, or Robin".to_string(),
                });
            }
        };
        let value = uniform_value(row)?;
        if value.dims != expected {
            return Err(ConductionError::Dimensions {
                context: format!("scenario row on {:?} ({:?})", row.region, row.kind),
                expected: expected.0,
                found: value.dims.0,
            });
        }
        match row.kind {
            BcKind::Dirichlet => {
                debug_assert_eq!(expected, TEMPERATURE_DIMS);
                ThermalBc::dirichlet(value.value)
            }
            BcKind::Neumann => {
                debug_assert_eq!(expected, HEAT_FLUX_DIMS);
                ThermalBc::neumann(value.value)
            }
            BcKind::Robin => {
                debug_assert_eq!(expected, HTC_DIMS);
                let reference = match t_ref {
                    Some(t) => t,
                    None => ambient_reference(row, environment)?,
                };
                ThermalBc::robin(value.value, reference)
            }
            other => Err(ConductionError::ScenarioRow {
                region: row.region.clone(),
                what: format!("boundary kind {other:?} is not a conduction row"),
                fix: "use Dirichlet, Neumann, or Robin".to_string(),
            }),
        }
    }
}

/// The uniform scalar an `fs-scenario` thermal row carries, or a typed
/// refusal naming the carrier this steady solve cannot consume.
fn uniform_value(row: &BoundaryCondition) -> Result<fs_qty::QtyAny, ConductionError> {
    match &row.value {
        Some(BcValue::Uniform(q)) => Ok(*q),
        Some(BcValue::Signal(_)) => Err(ConductionError::ScenarioRow {
            region: row.region.clone(),
            what: "the row carries a time signal".to_string(),
            fix: "this crate is STEADY-only; evaluate the signal at the instant being \
                  solved and pass a uniform value"
                .to_string(),
        }),
        Some(BcValue::Profile(_)) => Err(ConductionError::ScenarioRow {
            region: row.region.clone(),
            what: "the row carries a Chebyshev spatial profile".to_string(),
            fix: "sample the profile at the mesh's boundary vertices and pass a \
                  ScalarField::Nodal"
                .to_string(),
        }),
        Some(BcValue::Typed(payload)) => Err(ConductionError::ScenarioRow {
            region: row.region.clone(),
            what: format!("the row carries a typed {:?} payload", payload.kind()),
            fix: "thermal rows use the legacy scalar carrier in fs-scenario's own \
                  expectation table"
                .to_string(),
        }),
        None => Err(ConductionError::ScenarioRow {
            region: row.region.clone(),
            what: "the row carries no value".to_string(),
            fix: "supply the prescribed quantity".to_string(),
        }),
    }
}

/// The scenario ambient temperature, dimension-checked, used as a Robin
/// row's `T_ref` when the caller names no other reference.
fn ambient_reference(
    row: &BoundaryCondition,
    environment: &Environment,
) -> Result<f64, ConductionError> {
    let ambient = environment.ambient_temperature;
    if ambient.dims == TEMPERATURE_DIMS {
        Ok(ambient.value)
    } else {
        Err(ConductionError::Dimensions {
            context: format!(
                "scenario ambient temperature used as T_ref for {:?}",
                row.region
            ),
            expected: TEMPERATURE_DIMS.0,
            found: ambient.dims.0,
        })
    }
}

/// A complete boundary partition with one condition per region.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalBoundary {
    names: Vec<String>,
    conditions: Vec<ThermalBc>,
    /// Region index per boundary-face slot; `None` only when an explicit
    /// adiabatic remainder was declared.
    face_region: Vec<Option<usize>>,
    /// Prescribed vertices, ascending by vertex id.
    dirichlet: Vec<(usize, f64)>,
    adiabatic_remainder: usize,
}

impl ThermalBoundary {
    /// Region names in declaration order.
    #[must_use]
    pub fn region_names(&self) -> &[String] {
        &self.names
    }

    /// Conditions in declaration order.
    #[must_use]
    pub fn conditions(&self) -> &[ThermalBc] {
        &self.conditions
    }

    /// The condition owning a boundary-face slot, if any.
    #[must_use]
    pub fn condition_for(&self, boundary_face_slot: usize) -> Option<&ThermalBc> {
        self.face_region[boundary_face_slot].map(|r| &self.conditions[r])
    }

    /// The prescribed `(vertex, temperature)` pairs, ascending by vertex.
    #[must_use]
    pub fn dirichlet(&self) -> &[(usize, f64)] {
        &self.dirichlet
    }

    /// How many boundary faces fell to the declared adiabatic remainder.
    #[must_use]
    pub const fn adiabatic_remainder_faces(&self) -> usize {
        self.adiabatic_remainder
    }

    /// True when at least one region is a Robin row.
    #[must_use]
    pub fn has_robin(&self) -> bool {
        self.conditions
            .iter()
            .any(|c| matches!(c, ThermalBc::Robin { .. }))
    }
}

/// Builds a [`ThermalBoundary`] by tagging boundary faces with named
/// regions. Regions must PARTITION the boundary: a face claimed twice is
/// a refusal, and leftover faces are a refusal unless the caller says
/// [`ThermalBoundaryBuilder::adiabatic_remainder`] out loud.
pub struct ThermalBoundaryBuilder<'m> {
    mesh: &'m ConductionMesh,
    names: Vec<String>,
    conditions: Vec<ThermalBc>,
    face_region: Vec<Option<usize>>,
    adiabatic_remainder: bool,
}

impl<'m> ThermalBoundaryBuilder<'m> {
    /// Start from a fully untagged boundary.
    #[must_use]
    pub fn new(mesh: &'m ConductionMesh) -> ThermalBoundaryBuilder<'m> {
        ThermalBoundaryBuilder {
            mesh,
            names: Vec::new(),
            conditions: Vec::new(),
            face_region: vec![None; mesh.boundary().len()],
            adiabatic_remainder: false,
        }
    }

    /// Declare a region: every UNTAGGED boundary face satisfying
    /// `select` joins it and carries `condition`.
    ///
    /// # Errors
    /// [`ConductionError::DuplicateRegion`] for a repeated name;
    /// [`ConductionError::OverlappingRegion`] when `select` matches a
    /// face another region already owns; the condition's own validation
    /// refusals.
    pub fn region(
        mut self,
        name: &str,
        select: impl Fn(&BoundaryFace) -> bool,
        condition: ThermalBc,
    ) -> Result<ThermalBoundaryBuilder<'m>, ConductionError> {
        if self.names.iter().any(|n| n == name) {
            return Err(ConductionError::DuplicateRegion {
                region: name.to_string(),
            });
        }
        condition.validate(self.mesh.vertex_count())?;
        let index = self.names.len();
        for (slot, face) in self.mesh.boundary().iter().enumerate() {
            if !select(face) {
                continue;
            }
            if let Some(owner) = self.face_region[slot] {
                return Err(ConductionError::OverlappingRegion {
                    region: name.to_string(),
                    owner: self.names[owner].clone(),
                    face: slot,
                });
            }
            self.face_region[slot] = Some(index);
        }
        self.names.push(name.to_string());
        self.conditions.push(condition);
        Ok(self)
    }

    /// Declare a region holding every face still untagged. This is the
    /// explicit "everything else" clause: [`ThermalBoundaryBuilder::region`]
    /// refuses to steal a face another region already owns, so a
    /// catch-all predicate cannot silently override an earlier row.
    ///
    /// # Errors
    /// [`ConductionError::DuplicateRegion`]; the condition's own
    /// validation refusals.
    pub fn remainder(
        mut self,
        name: &str,
        condition: ThermalBc,
    ) -> Result<ThermalBoundaryBuilder<'m>, ConductionError> {
        if self.names.iter().any(|n| n == name) {
            return Err(ConductionError::DuplicateRegion {
                region: name.to_string(),
            });
        }
        condition.validate(self.mesh.vertex_count())?;
        let index = self.names.len();
        for slot in &mut self.face_region {
            if slot.is_none() {
                *slot = Some(index);
            }
        }
        self.names.push(name.to_string());
        self.conditions.push(condition);
        Ok(self)
    }

    /// Declare that every remaining untagged boundary face is adiabatic.
    /// This is the explicit form of the natural boundary condition — the
    /// model states it rather than inheriting it from silence.
    #[must_use]
    pub fn adiabatic_remainder(mut self) -> ThermalBoundaryBuilder<'m> {
        self.adiabatic_remainder = true;
        self
    }

    /// Finish the partition.
    ///
    /// # Errors
    /// [`ConductionError::UntaggedBoundary`] when faces are left over
    /// and no adiabatic remainder was declared.
    pub fn finish(self) -> Result<ThermalBoundary, ConductionError> {
        let untagged: Vec<usize> = self
            .face_region
            .iter()
            .enumerate()
            .filter_map(|(slot, r)| r.is_none().then_some(slot))
            .collect();
        if !untagged.is_empty() && !self.adiabatic_remainder {
            return Err(ConductionError::UntaggedBoundary {
                count: untagged.len(),
                first: untagged[0],
            });
        }

        // Dirichlet vertices: the union of the vertices of every face in
        // a Dirichlet region. Deterministic tie-break at shared vertices
        // (edges between two Dirichlet regions): the LOWEST declared
        // region index wins, so the pinned value never depends on face
        // traversal order.
        let mut assigned: Vec<Option<(usize, f64)>> = vec![None; self.mesh.vertex_count()];
        for (slot, region) in self.face_region.iter().enumerate() {
            let Some(region) = *region else { continue };
            let ThermalBc::Dirichlet { temperature } = &self.conditions[region] else {
                continue;
            };
            for &v in &self.mesh.boundary()[slot].vertices {
                let v = v as usize;
                let value = temperature.at(v);
                match assigned[v] {
                    None => assigned[v] = Some((region, value)),
                    Some((owner, _)) if region < owner => assigned[v] = Some((region, value)),
                    Some(_) => {}
                }
            }
        }
        let dirichlet: Vec<(usize, f64)> = assigned
            .iter()
            .enumerate()
            .filter_map(|(v, a)| a.map(|(_, value)| (v, value)))
            .collect();

        Ok(ThermalBoundary {
            names: self.names,
            conditions: self.conditions,
            face_region: self.face_region,
            dirichlet,
            adiabatic_remainder: untagged.len(),
        })
    }
}
