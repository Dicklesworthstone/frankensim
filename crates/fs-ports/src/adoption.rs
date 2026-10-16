//! Common module adoption and lossless N-port relation descriptions.
//!
//! These are schema contracts, not integrators. In particular, a skew matrix
//! establishes only interface power cancellation, never component passivity.

use crate::{
    ConservativeJunction, CoupleError, IntegratedWindowTransfer, PortSchema, StableId,
    WindowAuditReport,
};

/// Conversion between an owner's coordinates and its canonical SI port schema.
#[derive(Clone, Debug, PartialEq)]
pub enum PortCoordinateChart {
    /// Values already use the schema's coherent SI effort and flow.
    CanonicalSi,
    /// Mass-normalized modal coordinates, NOT point force/velocity samples.
    /// For reference mass M, F = f_q sqrt(M), v = q_dot / sqrt(M).
    /// Raw owner units remain N/sqrt(kg) and m sqrt(kg)/s, respectively.
    MassNormalized {
        /// Explicit coordinate reference mass [kg], not a specimen modal mass.
        reference_mass_kg: f64,
        /// Optional original signed swept-volume projection [m²/sqrt(kg)].
        /// This labels the pressure-work pullback; it is not a point-force area.
        volume_projection_m2_per_sqrt_kg: Option<f64>,
    },
}

/// Whether a module can supply a physical accounting-window energy audit.
#[derive(Clone, Debug, PartialEq)]
pub enum ModuleWindowAudit {
    /// The module has no accepted physical step yet.
    NotStepped,
    /// A complete retained account; violations remain visible in the report.
    Available(Box<WindowAuditReport>),
    /// A named owner or physical acceptance input is missing.
    Unavailable {
        /// Module whose accepted step is being described.
        owner: StableId,
        /// Specific missing owner, tolerance, or refused transfer input.
        reason: String,
    },
}

/// A module's ordered schemas and actual, accepted physical-time transfers.
pub trait PortedModule {
    /// Stable order: index i identifies the same port in every accepted window.
    fn port_schemas(&self) -> &[PortSchema];

    /// One explicit coordinate chart per schema, in the same stable order.
    fn port_coordinate_charts(&self) -> &[PortCoordinateChart];

    /// Internal N-way relations, after the owner's declared power charts.
    fn junction_schemas(&self) -> &[JunctionSchema];

    /// Physical duration used to pre-integrate the last accepted exchanges [s].
    fn last_window_duration_s(&self) -> Option<f64>;

    /// Last accepted step, integrated in physical time, positive into the owner.
    ///
    /// Held effort uses effort times the actual conjugate coordinate increment.
    /// Rejected trials must preserve this record. Before any accepted step this
    /// returns `None`; it must never substitute a fabricated zero-step account.
    fn last_window_transfers(&self) -> Option<&[IntegratedWindowTransfer]>;

    /// A failed observation never fabricates a zero transfer or changes motion.
    fn last_window_transfer_error(&self) -> Option<&CoupleError>;

    /// Audit the accepted account, or name the unavailable physical input.
    fn audit_last_window(&self) -> Result<ModuleWindowAudit, CoupleError>;
}

/// Admission ceiling, including every endpoint of a relation.
pub const MAX_JUNCTION_PORTS: usize = 4097;

/// A lossless relation in the schemas' common, power-conjugate coordinates.
#[derive(Clone, Debug, PartialEq)]
pub enum JunctionLaw {
    /// `e_i = s_i e`, `sum(s_i f_i) = 0`, with each sign equal to ±1.
    /// The `[1, 1]` two-port instance is the original conservative junction.
    CommonEffort {
        /// Effort pullback and matching flow-conservation signs.
        flow_signs: Vec<i8>,
    },
    /// `f_i = s_i f`, `sum(s_i e_i) = 0`, with each sign equal to ±1.
    CommonFlow {
        /// Flow pullback and matching effort-conservation signs.
        effort_signs: Vec<i8>,
    },
    /// `f = J e` in declared SI coordinates, with exactly skew `J`.
    /// Coefficients carry the constitutive flow/effort conversion of these
    /// common coordinates; this does not supply a cross-domain area chart.
    Dirac {
        /// Square matrix order, equal to the number of endpoints.
        rows: usize,
        /// Row-major J. Diagonal entries are zero and J_ij = -J_ji exactly.
        incidence: Vec<f64>,
    },
}

/// One N-way lossless relation, with identities and port ordering retained.
#[derive(Clone, Debug, PartialEq)]
pub struct JunctionSchema {
    id: StableId,
    ports: Vec<PortSchema>,
    law: JunctionLaw,
}

impl JunctionSchema {
    /// Admit conjugate endpoints and an algebraically power-cancelling law.
    ///
    /// Endpoint compatibility uses exactly [`ConservativeJunction::new`],
    /// including coordinates, orientation, clocks and conservation roles.
    /// Mixed physical coordinates require an explicit upstream power chart;
    /// they cannot be made conjugate by relabelling their dimensions.
    ///
    /// # Errors
    /// Refuses fewer than two or more than [`MAX_JUNCTION_PORTS`] endpoints,
    /// identity aliases, incompatible schemas, malformed sign vectors, and
    /// nonfinite or non-skew matrices. Structural equality has zero tolerance.
    pub fn try_new(
        id: StableId,
        ports: Vec<PortSchema>,
        law: JunctionLaw,
    ) -> Result<Self, CoupleError> {
        let n = ports.len();
        if !(2..=MAX_JUNCTION_PORTS).contains(&n) {
            return Err(CoupleError::InvalidJunction {
                field: "port_count",
            });
        }
        for i in 1..n {
            ConservativeJunction::new(id.clone(), ports[0].clone(), ports[i].clone())?;
            for earlier in &ports[1..i] {
                if earlier.id() == ports[i].id() {
                    return Err(CoupleError::DuplicateIdentity {
                        id: ports[i].id().as_str().to_owned(),
                    });
                }
            }
        }
        match &law {
            JunctionLaw::CommonEffort { flow_signs: signs }
            | JunctionLaw::CommonFlow {
                effort_signs: signs,
            } => {
                if signs.len() != n || signs.iter().any(|&s| s != 1 && s != -1) {
                    return Err(CoupleError::InvalidJunction { field: "signs" });
                }
            }
            JunctionLaw::Dirac { rows, incidence } => {
                if *rows != n || incidence.len() != n * n {
                    return Err(CoupleError::InvalidJunction {
                        field: "matrix_shape",
                    });
                }
                for i in 0..n {
                    if incidence[i * n + i] != 0.0 {
                        return Err(CoupleError::InvalidJunction {
                            field: "skew_diagonal",
                        });
                    }
                    for j in 0..n {
                        let a = incidence[i * n + j];
                        let b = incidence[j * n + i];
                        if !a.is_finite() || !b.is_finite() || a != -b {
                            return Err(CoupleError::InvalidJunction {
                                field: "skew_incidence",
                            });
                        }
                    }
                }
            }
        }
        Ok(Self { id, ports, law })
    }

    /// Caller-supplied relation identity.
    pub fn id(&self) -> &StableId {
        &self.id
    }
    /// Endpoints in original module/graph order.
    pub fn ports(&self) -> &[PortSchema] {
        &self.ports
    }
    /// Admitted relation in those coordinates.
    pub fn law(&self) -> &JunctionLaw {
        &self.law
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp};

    fn id(s: &str) -> StableId {
        StableId::new(s).unwrap()
    }
    fn ports(n: usize) -> Vec<PortSchema> {
        (0..n)
            .map(|i| {
                PortKind::MechanicalForceVelocity
                    .scalar_seed_schema(
                        id(&format!("body/{i}")),
                        CoordinateBinding::new(
                            id("cartesian"),
                            id("wall"),
                            PortOrientation::OutwardFromOwner,
                        ),
                        PortTimestamp::new(id("mechanical-clock"), 0),
                    )
                    .unwrap()
            })
            .collect()
    }
    fn effort(n: usize) -> JunctionLaw {
        JunctionLaw::CommonEffort {
            flow_signs: vec![1; n],
        }
    }

    #[test]
    fn g0_two_three_sixty_four_and_cap_are_admitted() {
        for n in [2, 3, 64, MAX_JUNCTION_PORTS] {
            let relation = JunctionSchema::try_new(id("junction"), ports(n), effort(n)).unwrap();
            assert_eq!(relation.ports().len(), n);
            JunctionSchema::try_new(
                id("junction"),
                ports(n),
                JunctionLaw::CommonFlow {
                    effort_signs: vec![1; n],
                },
            )
            .unwrap();
        }
        for n in [0, 1, MAX_JUNCTION_PORTS + 1] {
            assert_eq!(
                JunctionSchema::try_new(id("junction"), ports(n), effort(n)),
                Err(CoupleError::InvalidJunction {
                    field: "port_count"
                })
            );
        }
    }

    #[test]
    fn g0_two_port_admission_matches_original_including_errors() {
        let p = ports(2);
        let mut cases = vec![(id("junction"), p.clone()), (p[0].id().clone(), p.clone())];
        let mut duplicate = p.clone();
        duplicate[1] = duplicate[0].clone();
        cases.push((id("junction"), duplicate));
        for orientation in [PortOrientation::AlongFrame, PortOrientation::AgainstFrame] {
            let mut q = p.clone();
            q[1].coordinates.orientation = orientation;
            cases.push((id("junction"), q));
        }
        let mut q = p.clone();
        q[1].timestamp = PortTimestamp::new(id("other-clock"), 0);
        cases.push((id("junction"), q));
        let mut q = p.clone();
        q[1].coordinates.basis = id("other-basis");
        cases.push((id("junction"), q));
        let mut q = p.clone();
        q[1].coordinates.frame = id("other-frame");
        cases.push((id("junction"), q));
        let mut q = p.clone();
        q[1] = PortKind::FluidPressureFlux
            .scalar_seed_schema(
                id("fluid"),
                p[0].coordinates().clone(),
                p[0].timestamp().clone(),
            )
            .unwrap();
        cases.push((id("junction"), q));
        let mut q = p.clone();
        q[1].conservation_roles
            .push(crate::ConservationRole::LinearMomentum);
        cases.push((id("junction"), q));
        for (name, p) in cases {
            let old =
                ConservativeJunction::new(name.clone(), p[0].clone(), p[1].clone()).map(|_| ());
            let new = JunctionSchema::try_new(name, p, effort(2)).map(|_| ());
            assert_eq!(new, old);
        }
    }

    #[test]
    fn g0_skew_relation_cancels_power_and_rejects_mutations() {
        let matrix = vec![0.0, 2.0, -3.0, -2.0, 0.0, 4.0, 3.0, -4.0, 0.0];
        JunctionSchema::try_new(
            id("junction"),
            ports(3),
            JunctionLaw::Dirac {
                rows: 3,
                incidence: matrix.clone(),
            },
        )
        .unwrap();
        let e = [5.0, 7.0, 11.0];
        let power: f64 = (0..3)
            .map(|i| e[i] * (0..3).map(|j| matrix[3 * i + j] * e[j]).sum::<f64>())
            .sum();
        assert_eq!(power, 0.0);
        for (index, value) in [(0, 1.0), (1, -2.0), (1, f64::NAN), (1, f64::INFINITY)] {
            let mut bad = matrix.clone();
            bad[index] = value;
            assert!(
                JunctionSchema::try_new(
                    id("junction"),
                    ports(3),
                    JunctionLaw::Dirac {
                        rows: 3,
                        incidence: bad
                    }
                )
                .is_err()
            );
        }
        for law in [
            JunctionLaw::CommonEffort {
                flow_signs: vec![1, 0, 1],
            },
            JunctionLaw::CommonFlow {
                effort_signs: vec![1, 1],
            },
            JunctionLaw::Dirac {
                rows: 2,
                incidence: matrix,
            },
        ] {
            assert!(JunctionSchema::try_new(id("junction"), ports(3), law).is_err());
        }
        let mut p = ports(3);
        p[2] = p[1].clone();
        assert!(JunctionSchema::try_new(id("junction"), p, effort(3)).is_err());
    }
}
