//! kagp increment-1 conformance: law-node admission, the fs-matdb
//! registry, consistent-tangent gates, and the aggregate state codec.
//! Every rejection logs node id, law id, and the failed obligation.

use std::collections::BTreeMap;

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_material::graph::{
    AggregateStateSchema, Differentiability, EnergyBehavior, GraphError, LawNode, LawRegistry,
    NodeDeclaration, NodeOutput, NodeRole, Port, TimeParity, admit_node, check_consistent_tangent,
    check_free_energy_consistency, check_onsager_casimir, check_psd_symmetric_part,
};
use fs_qty::Dims;

const CONDUCTIVITY_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);
const GRAD_T_DIMS: Dims = Dims([-1, 0, 0, 1, 0, 0]);
const FLUX_DIMS: Dims = Dims([0, 1, -3, 0, 0, 0]);

fn fourier_card(k: f64) -> ConstitutiveModelCard {
    let mut parameters = BTreeMap::new();
    parameters.insert(
        "conductivity".to_string(),
        LawParameter {
            value: k,
            dims: CONDUCTIVITY_DIMS,
        },
    );
    ConstitutiveModelCard {
        law: LawId("fourier-conduction".to_string()),
        law_version: 1,
        parameters,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
        sources: vec![hash_bytes(b"conductivity calibration")],
        provenance: Provenance {
            source: "handbook".to_string(),
            license: "internal-use".to_string(),
            artifact: None,
        },
    }
}

/// Fixture bulk-transport node: q = -k * gradT (1-D Fourier), with a
/// correct analytic tangent.
struct FourierNode {
    declaration: NodeDeclaration,
    k: f64,
}

impl FourierNode {
    fn from_card(card: &ConstitutiveModelCard) -> Result<Box<dyn LawNode>, GraphError> {
        let k = card
            .parameters
            .get("conductivity")
            .ok_or_else(|| GraphError::MissingParameter {
                law: card.law.0.clone(),
                parameter: "conductivity".to_string(),
            })?
            .value;
        Ok(Box::new(FourierNode {
            declaration: NodeDeclaration {
                law: card.law.clone(),
                law_version: card.law_version,
                role: NodeRole::BulkTransport,
                inputs: vec![Port {
                    name: "grad_T".to_string(),
                    dims: GRAD_T_DIMS,
                    parity: TimeParity::Even,
                }],
                outputs: vec![Port {
                    name: "heat_flux".to_string(),
                    dims: FLUX_DIMS,
                    parity: TimeParity::Even,
                }],
                state_slots: Vec::new(),
                state_schema_version: card.state_schema_version,
                calibration: card.validity.clone(),
                differentiability: Differentiability::Smooth,
                energy: EnergyBehavior::NonNegativeDissipation,
                tangent_claimed: true,
            },
            k,
        }))
    }
}

impl LawNode for FourierNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }

    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        let grad_t = inputs[0];
        let flux = -self.k * grad_t;
        Ok(NodeOutput {
            outputs: vec![flux],
            next_state: Vec::new(),
            // sigma = q * (-gradT) / T ... reported here as k*gradT^2
            // scaled: non-negative for k >= 0 by construction.
            dissipation_rate: Some(self.k * grad_t * grad_t),
        })
    }

    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![-self.k])
    }
}

/// A node that CLAIMS a tangent but lies about it (wrong sign).
struct WrongTangentNode(FourierNode);

impl LawNode for WrongTangentNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.0.declaration
    }
    fn evaluate(&self, state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        self.0.evaluate(state, inputs)
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![self.0.k])
    }
}

fn built_fourier() -> Box<dyn LawNode> {
    FourierNode::from_card(&fourier_card(40.0)).expect("fixture builds")
}

#[test]
fn admission_rejects_incomplete_declarations_with_typed_diagnostics() {
    let good = built_fourier();
    admit_node("bulk/steel", good.as_ref()).expect("complete declaration admits");

    struct Broken {
        declaration: NodeDeclaration,
    }
    impl LawNode for Broken {
        fn declaration(&self) -> &NodeDeclaration {
            &self.declaration
        }
        fn evaluate(&self, _s: &[f64], _i: &[f64]) -> Result<NodeOutput, GraphError> {
            Ok(NodeOutput {
                outputs: vec![0.0],
                next_state: Vec::new(),
                dissipation_rate: None,
            })
        }
    }

    let template = built_fourier().declaration().clone();

    let no_outputs = Broken {
        declaration: NodeDeclaration {
            outputs: Vec::new(),
            ..template.clone()
        },
    };
    let refusal = admit_node("bulk/broken", &no_outputs).expect_err("no outputs refuses");
    assert!(matches!(
        &refusal,
        GraphError::IncompleteDeclaration { node, law, .. }
            if node == "bulk/broken" && law == "fourier-conduction"
    ));
    println!("{{\"case\":\"admission-refusal\",\"log\":\"{refusal}\"}}");

    let dup_ports = Broken {
        declaration: NodeDeclaration {
            inputs: vec![
                Port {
                    name: "grad_T".to_string(),
                    dims: GRAD_T_DIMS,
                    parity: TimeParity::Even,
                },
                Port {
                    name: "grad_T".to_string(),
                    dims: GRAD_T_DIMS,
                    parity: TimeParity::Even,
                },
            ],
            ..template.clone()
        },
    };
    assert!(matches!(
        admit_node("bulk/dup", &dup_ports),
        Err(GraphError::IncompleteDeclaration { .. })
    ));

    let claims_unmet = Broken {
        declaration: NodeDeclaration {
            tangent_claimed: true,
            ..template.clone()
        },
    };
    assert!(matches!(
        admit_node("bulk/liar", &claims_unmet),
        Err(GraphError::TangentClaimUnmet { .. })
    ));

    let external = Broken {
        declaration: NodeDeclaration {
            role: NodeRole::TopologyBalance,
            ..template
        },
    };
    assert!(matches!(
        admit_node("bulk/external", &external),
        Err(GraphError::ExternallyOwnedRole { .. })
    ));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"node-admission\",\"verdict\":\"pass\",\
         \"detail\":\"missing outputs, duplicate ports, unmet tangent claims, and external roles \
         refuse with node+law+obligation\"}}"
    );
}

#[test]
fn registry_instantiates_from_validated_cards_and_refuses_drift() {
    let mut registry = LawRegistry::new();
    registry.register(&LawId("fourier-conduction".to_string()), 1, |card| {
        FourierNode::from_card(card)
    });

    let node = registry
        .instantiate("bulk/steel", &fourier_card(40.0))
        .expect("registered card instantiates");
    let out = node.evaluate(&[], &[5.0]).expect("evaluates");
    assert_eq!(out.outputs, vec![-200.0]);

    let mut unknown_version = fourier_card(40.0);
    unknown_version.law_version = 2;
    assert!(matches!(
        registry.instantiate("bulk/steel", &unknown_version),
        Err(GraphError::UnknownLaw { version: 2, .. })
    ));

    let mut missing_parameter = fourier_card(40.0);
    missing_parameter.parameters.clear();
    missing_parameter.parameters.insert(
        "not-conductivity".to_string(),
        LawParameter {
            value: 1.0,
            dims: CONDUCTIVITY_DIMS,
        },
    );
    assert!(matches!(
        registry.instantiate("bulk/steel", &missing_parameter),
        Err(GraphError::MissingParameter { .. })
    ));

    let mut invalid_card = fourier_card(f64::NAN);
    invalid_card.law_version = 1;
    assert!(matches!(
        registry.instantiate("bulk/steel", &invalid_card),
        Err(GraphError::Card(_))
    ));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"registry\",\"verdict\":\"pass\",\
         \"detail\":\"cards instantiate through validation; version/parameter/card drift refuses\"}}"
    );
}

#[test]
fn consistent_tangent_gate_passes_honest_nodes_and_catches_liars() {
    let honest = built_fourier();
    check_consistent_tangent("bulk/steel", honest.as_ref(), &[], &[3.0], 1e-6)
        .expect("honest tangent passes the FD gate");

    let liar = WrongTangentNode(FourierNode {
        declaration: built_fourier().declaration().clone(),
        k: 40.0,
    });
    let refusal = check_consistent_tangent("bulk/liar", &liar, &[], &[3.0], 1e-6)
        .expect_err("wrong-sign tangent fails the gate");
    assert!(matches!(
        &refusal,
        GraphError::IncompleteDeclaration { obligation, .. }
            if obligation.contains("finite-difference")
    ));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"tangent-gate\",\"verdict\":\"pass\",\
         \"detail\":\"FD probe accepts the honest tangent and names the liar's divergence\"}}"
    );
}

#[test]
fn aggregate_state_codec_round_trips_and_refuses_drift() {
    let stateless = built_fourier().declaration().clone();
    let mut memory = built_fourier().declaration().clone();
    memory.state_slots = vec!["eps_p".to_string(), "alpha".to_string()];
    memory.state_schema_version = 3;

    let schema =
        AggregateStateSchema::assemble(&[("bulk/steel", &stateless), ("memory/voce", &memory)]);
    assert_eq!(schema.total_slots(), 2);

    let (version, buffer) = schema
        .encode(&[&[], &[0.01, 150.0e6]])
        .expect("encode succeeds");
    let decoded = schema.decode(version, &buffer).expect("round trip");
    assert_eq!(decoded, vec![Vec::new(), vec![0.01, 150.0e6]]);

    assert!(matches!(
        schema.decode(version ^ 1, &buffer),
        Err(GraphError::StateSchemaMismatch {
            obligation: "schema version differs",
            ..
        })
    ));
    assert!(matches!(
        schema.decode(version, &buffer[..1]),
        Err(GraphError::StateSchemaMismatch {
            obligation: "buffer length differs",
            ..
        })
    ));
    assert!(matches!(
        schema.encode(&[&[]]),
        Err(GraphError::StateSchemaMismatch {
            obligation: "node count differs",
            ..
        })
    ));

    // ANY layout change moves the schema version.
    let mut renamed = memory.clone();
    renamed.state_slots = vec!["eps_p".to_string(), "beta".to_string()];
    let other =
        AggregateStateSchema::assemble(&[("bulk/steel", &stateless), ("memory/voce", &renamed)]);
    assert_ne!(schema.version(), other.version());
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"state-codec\",\"verdict\":\"pass\",\
         \"detail\":\"round trip exact; version/length/count drift refuses; layout moves version\"}}"
    );
}

/// Increment-2 fixture: a linear-elastic storage node, psi = k/2 x².
/// `gradient_scale` != 1 makes the outputs LIE about being conjugate
/// forces (the energy-consistency gate must catch it).
struct ElasticEnergyNode {
    declaration: NodeDeclaration,
    k: f64,
    gradient_scale: f64,
    provide_energy: bool,
}

impl ElasticEnergyNode {
    fn new(k: f64, gradient_scale: f64, provide_energy: bool) -> ElasticEnergyNode {
        ElasticEnergyNode {
            declaration: NodeDeclaration {
                law: LawId("linear-elastic-energy".to_string()),
                law_version: 1,
                role: NodeRole::InternalMemory,
                inputs: vec![Port {
                    name: "strain".to_string(),
                    dims: Dims([0, 0, 0, 0, 0, 0]),
                    parity: TimeParity::Even,
                }],
                outputs: vec![Port {
                    name: "stress".to_string(),
                    dims: Dims([-1, 1, -2, 0, 0, 0]),
                    parity: TimeParity::Even,
                }],
                state_slots: Vec::new(),
                state_schema_version: 1,
                calibration: ValidityDomain::unconstrained(),
                differentiability: Differentiability::Smooth,
                energy: EnergyBehavior::FreeEnergyStorage,
                tangent_claimed: true,
            },
            k,
            gradient_scale,
            provide_energy,
        }
    }
}

impl LawNode for ElasticEnergyNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        Ok(NodeOutput {
            outputs: vec![self.gradient_scale * self.k * inputs[0]],
            next_state: Vec::new(),
            dissipation_rate: None,
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(vec![self.gradient_scale * self.k])
    }
    fn free_energy(&self, _state: &[f64], inputs: &[f64]) -> Option<f64> {
        self.provide_energy
            .then(|| 0.5 * self.k * inputs[0] * inputs[0])
    }
}

/// Increment-2 fixture: a 2x2 coupling block with declared parities and
/// a controllable off-diagonal relationship.
struct CoupledBlockNode {
    declaration: NodeDeclaration,
    matrix: [f64; 4],
}

impl CoupledBlockNode {
    fn new(matrix: [f64; 4], parities: [TimeParity; 2]) -> CoupledBlockNode {
        let port = |name: &str, parity: TimeParity| Port {
            name: name.to_string(),
            dims: Dims([0, 0, 0, 0, 0, 0]),
            parity,
        };
        CoupledBlockNode {
            declaration: NodeDeclaration {
                law: LawId("coupled-block".to_string()),
                law_version: 1,
                role: NodeRole::ReversibleCoupling,
                inputs: vec![port("force_a", parities[0]), port("force_b", parities[1])],
                outputs: vec![port("flux_a", parities[0]), port("flux_b", parities[1])],
                state_slots: Vec::new(),
                state_schema_version: 1,
                calibration: ValidityDomain::unconstrained(),
                differentiability: Differentiability::Smooth,
                energy: EnergyBehavior::NonNegativeDissipation,
                tangent_claimed: true,
            },
            matrix,
        }
    }
}

impl LawNode for CoupledBlockNode {
    fn declaration(&self) -> &NodeDeclaration {
        &self.declaration
    }
    fn evaluate(&self, _state: &[f64], inputs: &[f64]) -> Result<NodeOutput, GraphError> {
        Ok(NodeOutput {
            outputs: vec![
                self.matrix[0] * inputs[0] + self.matrix[1] * inputs[1],
                self.matrix[2] * inputs[0] + self.matrix[3] * inputs[1],
            ],
            next_state: Vec::new(),
            dissipation_rate: None,
        })
    }
    fn tangent(&self, _state: &[f64], _inputs: &[f64]) -> Option<Vec<f64>> {
        Some(self.matrix.to_vec())
    }
}

#[test]
fn free_energy_gates_pass_storage_laws_and_catch_violations() {
    let honest = ElasticEnergyNode::new(2.0e9, 1.0, true);
    admit_node("memory/elastic", &honest).expect("storage node with psi admits");
    check_free_energy_consistency("memory/elastic", &honest, &[], &[1e-3], 1e-3)
        .expect("conjugate-force gradient and symmetric Hessian pass");

    let liar = ElasticEnergyNode::new(2.0e9, 2.0, true);
    let refusal = check_free_energy_consistency("memory/liar", &liar, &[], &[1e-3], 1e-3)
        .expect_err("scaled outputs are not the conjugate forces");
    assert!(matches!(
        &refusal,
        GraphError::IncompleteDeclaration { obligation, .. }
            if obligation.contains("conjugate forces")
    ));

    let no_psi = ElasticEnergyNode::new(2.0e9, 1.0, false);
    assert!(matches!(
        admit_node("memory/no-psi", &no_psi),
        Err(GraphError::EnergyClaimUnmet { .. })
    ));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"free-energy-gates\",\"verdict\":\"pass\",\
         \"detail\":\"conjugate-force+Hessian gates pass honest storage, catch scaled gradients \
         and missing psi\"}}"
    );
}

#[test]
fn psd_gate_enforces_the_second_law_on_transport_blocks() {
    // Force -> flux convention: sigma / kappa matrices must have a PSD
    // symmetric part.
    let passive = CoupledBlockNode::new([4.0, 0.3, 0.3, 2.0], [TimeParity::Even, TimeParity::Even]);
    check_psd_symmetric_part("bulk/passive", &passive, &[], &[0.0, 0.0], 1e-12)
        .expect("passive conductivity passes");

    let active = CoupledBlockNode::new([4.0, 0.3, 0.3, -2.0], [TimeParity::Even, TimeParity::Even]);
    let refusal = check_psd_symmetric_part("bulk/active", &active, &[], &[0.0, 0.0], 1e-12)
        .expect_err("negative branch conductivity fails");
    assert!(matches!(
        &refusal,
        GraphError::IncompleteDeclaration { obligation, .. }
            if obligation.contains("positive semidefinite")
    ));
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"psd-gate\",\"verdict\":\"pass\",\
         \"detail\":\"Sylvester audit of the symmetric part passes passive blocks, refuses \
         active ones\"}}"
    );
}

#[test]
fn onsager_casimir_gate_separates_reciprocity_classes() {
    // Even-even coupling (thermoelectric): off-diagonals must MATCH.
    let onsager = CoupledBlockNode::new([4.0, 0.3, 0.3, 2.0], [TimeParity::Even, TimeParity::Even]);
    check_onsager_casimir("couple/thermoelectric", &onsager, &[], &[0.0, 0.0], 1e-12)
        .expect("symmetric even-even coupling passes");
    let broken = CoupledBlockNode::new([4.0, 0.3, -0.3, 2.0], [TimeParity::Even, TimeParity::Even]);
    assert!(
        check_onsager_casimir("couple/broken", &broken, &[], &[0.0, 0.0], 1e-12).is_err(),
        "antisymmetric even-even coupling violates Onsager"
    );

    // Mixed parity (gyroscopic/Hall class): off-diagonals must be
    // ANTISYMMETRIC (Casimir), and the symmetric version must refuse.
    let casimir = CoupledBlockNode::new([4.0, 0.3, -0.3, 2.0], [TimeParity::Even, TimeParity::Odd]);
    check_onsager_casimir("couple/gyroscopic", &casimir, &[], &[0.0, 0.0], 1e-12)
        .expect("antisymmetric mixed-parity coupling passes");
    let wrong_class =
        CoupledBlockNode::new([4.0, 0.3, 0.3, 2.0], [TimeParity::Even, TimeParity::Odd]);
    assert!(
        check_onsager_casimir("couple/wrong-class", &wrong_class, &[], &[0.0, 0.0], 1e-12).is_err(),
        "symmetric mixed-parity coupling violates Casimir"
    );
    println!(
        "{{\"suite\":\"fs-material\",\"case\":\"onsager-casimir\",\"verdict\":\"pass\",\
         \"detail\":\"even-even demands symmetry, mixed parity demands antisymmetry; both \
         violations refuse\"}}"
    );
}
