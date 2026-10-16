//! Thin observations of existing accepted steps. No solver is called here.
//!
//! Module-local IDs and clock ordinals are stable. A graph must namespace IDs
//! when embedding multiple instances; an ordinal never supplies a duration.
//! Internal connection work is retained at its port but excluded from the
//! external boundary of a whole-module account. Losses and prescribed sources
//! are separately disclosed, using the existing step owner's measured values.

use fs_ports::*;
use fs_qty::{ElectricCharge, Energy, Mass};
use std::sync::OnceLock;

pub(crate) fn id(value: &str) -> StableId {
    StableId::new(value).expect("adapter identifiers use canonical static segments")
}

/// Scalar schema in an explicit local chart. Runtime values are not admitted
/// here: existing owners remain responsible for their physical inputs.
pub(crate) fn schema(owner: &str, name: &str, kind: PortKind) -> PortSchema {
    kind.scalar_seed_schema(
        id(&format!("{owner}/{name}")),
        CoordinateBinding::new(
            id("adapter/si"),
            id("adapter/local"),
            PortOrientation::OutwardFromOwner,
        ),
        PortTimestamp::new(id(&format!("{owner}/clock")), 0),
    )
    .expect("static scalar schema has the kind's canonical SI dimensions")
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ModulePorts {
    owner: String,
    pub(crate) schemas: Vec<PortSchema>,
    pub(crate) charts: Vec<PortCoordinateChart>,
    pub(crate) junctions: Vec<JunctionSchema>,
    last: Option<StepAccount>,
    raw: Vec<PortTransfer>,
    integrated: OnceLock<Result<Vec<IntegratedWindowTransfer>, CoupleError>>,
}

/// Actual signed work/heat on one port. `None` denotes a mechanical exchange;
/// a thermal exchange carries the owner's actual discrete effort temperature.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct PortTransfer {
    pub(crate) energy_j: f64,
    pub(crate) temperature_k: Option<f64>,
    pub(crate) boundary: bool,
}
impl PortTransfer {
    pub(crate) fn work(energy_j: f64, boundary: bool) -> Self {
        Self {
            energy_j,
            temperature_k: None,
            boundary,
        }
    }
    pub(crate) fn heat(energy_j: f64, temperature_k: f64) -> Self {
        Self {
            energy_j,
            temperature_k: Some(temperature_k),
            boundary: true,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct StepAccount {
    pub(crate) tick: u64,
    pub(crate) duration_s: f64,
    pub(crate) storage_change_j: Option<f64>,
    pub(crate) loss_j: Option<f64>,
    pub(crate) included_work_j: f64,
    pub(crate) entropy_change_j_k: f64,
    pub(crate) tolerance_j: Option<f64>,
    pub(crate) unavailable: Option<&'static str>,
}

impl ModulePorts {
    pub(crate) fn new(owner: &str, ports: Vec<(&str, PortKind)>) -> Self {
        let schemas: Vec<_> = ports
            .into_iter()
            .map(|(name, kind)| schema(owner, name, kind))
            .collect();
        let charts = vec![PortCoordinateChart::CanonicalSi; schemas.len()];
        let raw = Vec::with_capacity(schemas.len());
        Self {
            owner: owner.into(),
            schemas,
            charts,
            junctions: Vec::new(),
            last: None,
            raw,
            integrated: OnceLock::new(),
        }
    }
    pub(crate) fn mechanical(owner: &str, count: usize) -> Self {
        let names: Vec<_> = (0..count).map(|i| format!("connection/{i}")).collect();
        Self::new(
            owner,
            names
                .iter()
                .map(|s| (s.as_str(), PortKind::MechanicalForceVelocity))
                .collect(),
        )
    }
    pub(crate) fn next_tick(&self) -> u64 {
        self.last.as_ref().map_or(1, |s| s.tick.saturating_add(1))
    }
    pub(crate) fn duration_s(&self) -> Option<f64> {
        self.last.as_ref().map(|s| s.duration_s)
    }
    fn interval(&self, tick: u64) -> Result<AccountingWindowInterval, CoupleError> {
        AccountingWindowInterval::try_new(
            id(&format!("{}/window/{tick}", self.owner)),
            PortTimestamp::new(id(&format!("{}/clock", self.owner)), tick.saturating_sub(1)),
            PortTimestamp::new(id(&format!("{}/clock", self.owner)), tick),
        )
    }
    /// Publish pre-integrated scalar observations without allocating. Expanded
    /// evidence and typed transfers are built lazily when an observer asks.
    pub(crate) fn accept(
        &mut self,
        account: StepAccount,
        ports: impl IntoIterator<Item = PortTransfer>,
    ) {
        self.raw.clear();
        self.raw.extend(ports);
        self.last = Some(account);
        self.integrated = OnceLock::new();
    }
    /// An independent observation may refuse after the physical owner accepts.
    /// Retain that refusal without manufacturing transfers or changing motion.
    pub(crate) fn accept_result(
        &mut self,
        account: StepAccount,
        ports: Result<Vec<PortTransfer>, CoupleError>,
    ) {
        match ports {
            Ok(ports) => self.accept(account, ports),
            Err(error) => {
                self.accept(account, []);
                self.integrated = OnceLock::from(Err(error));
            }
        }
    }
    fn integrated(&self) -> Option<&Result<Vec<IntegratedWindowTransfer>, CoupleError>> {
        let account = self.last.as_ref()?;
        Some(self.integrated.get_or_init(|| {
            let interval = self.interval(account.tick)?;
            if self.raw.len() != self.schemas.len() {
                return Err(CoupleError::AccountingWindowContextMismatch {
                    field: "adapter_port_count",
                });
            }
            self.raw
                .iter()
                .enumerate()
                .map(|(i, p)| {
                    transfer(p.energy_j, p.temperature_k, &interval, &format!("port/{i}"))
                })
                .collect()
        }))
    }
    pub(crate) fn transfers(&self) -> Option<&[IntegratedWindowTransfer]> {
        self.integrated()?.as_ref().ok().map(|s| s.as_slice())
    }
    pub(crate) fn transfer_error(&self) -> Option<&CoupleError> {
        self.integrated()?.as_ref().err()
    }
    pub(crate) fn audit(&self) -> Result<ModuleWindowAudit, CoupleError> {
        let Some(last) = &self.last else {
            return Ok(ModuleWindowAudit::NotStepped);
        };
        if let Some(reason) = last.unavailable {
            return Ok(ModuleWindowAudit::Unavailable {
                owner: id(&self.owner),
                reason: reason.into(),
            });
        }
        let Some(tolerance) = last.tolerance_j else {
            return Ok(ModuleWindowAudit::Unavailable {
                owner: id(&self.owner),
                reason: "owner has no declared energy tolerance [J]".into(),
            });
        };
        let (Some(storage_change), Some(loss)) = (last.storage_change_j, last.loss_j) else {
            return Ok(ModuleWindowAudit::Unavailable {
                owner: id(&self.owner),
                reason: "owner has no measured storage change/loss [J]".into(),
            });
        };
        let interval = self.interval(last.tick)?;
        let transfers = self
            .integrated()
            .expect("accepted account")
            .as_ref()
            .map_err(Clone::clone)?;
        let reference = id("accepted-step/initial-energy-zero");
        let convention = id("accepted-step/midpoint-heat-entropy");
        let mut manifest = Vec::new();
        let mut values = Vec::new();
        for (i, original) in self
            .schemas
            .iter()
            .enumerate()
            .filter(|(i, _)| self.raw[*i].boundary)
        {
            let mut port = schema(&self.owner, &format!("audit/{i}"), original.kind());
            // PortSchema has no public timestamp setter: admit the new window binding.
            port = PortSchema::try_new(
                port.id().clone(),
                port.kind(),
                port.effort_dimensions(),
                port.flow_dimensions(),
                port.shape(),
                port.coordinates().clone(),
                port.power_pairing(),
                interval.end().clone(),
                port.conservation_roles().iter().copied(),
            )?;
            let key = format!("port/{i}");
            manifest.push(boundary_entry(&port, &key, &interval)?);
            values.push((key, transfers[i].clone()));
        }
        // Named external sinks and prescribed work remain measured owner values.
        for (key, energy) in [
            ("dissipation", -loss),
            ("prescribed-work", last.included_work_j),
        ] {
            let base = schema(&self.owner, key, PortKind::MechanicalForceVelocity);
            let port = PortSchema::try_new(
                base.id().clone(),
                base.kind(),
                base.effort_dimensions(),
                base.flow_dimensions(),
                base.shape(),
                base.coordinates().clone(),
                base.power_pairing(),
                interval.end().clone(),
                base.conservation_roles().iter().copied(),
            )?;
            manifest.push(boundary_entry(&port, key, &interval)?);
            values.push((key.into(), transfer(energy, None, &interval, key)?));
        }
        let charge = WindowChargeSchema::ProvenNeutral {
            projection_evidence: evidence(
                &interval,
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::ProvenNeutral,
                },
                "neutral",
            ),
        };
        let spec = AccountingWindowSpec::try_new(
            interval.clone(),
            reference.clone(),
            WindowElementSchema::not_applicable(id("closed-module/no-material-transfer")),
            charge,
            convention,
            manifest.clone(),
            WindowAuditTolerances::try_new(
                Energy::new(tolerance),
                Mass::new(0.0),
                fs_qty::Amount::new(0.0),
                ElectricCharge::new(0.0),
                Entropy::new(0.0),
            )?,
            None,
            evidence(
                &interval,
                WindowEvidenceRole::manifest_closure(&manifest),
                "closure",
            ),
        )?;
        let inventory = |endpoint, energy, entropy| {
            WindowInventorySnapshot::try_new(
                &spec,
                endpoint,
                WindowBalance::try_new(
                    Energy::new(energy),
                    Mass::new(0.0),
                    [],
                    ElectricCharge::new(0.0),
                    Entropy::new(entropy),
                )?,
                evidence(
                    &interval,
                    WindowEvidenceRole::Inventory {
                        endpoint,
                        energy_reference: reference.clone(),
                    },
                    if endpoint == WindowEndpoint::Initial {
                        "initial"
                    } else {
                        "final"
                    },
                ),
            )
        };
        let initial = inventory(WindowEndpoint::Initial, 0.0, 0.0)?;
        let final_snapshot = inventory(
            WindowEndpoint::Final,
            storage_change,
            last.entropy_change_j_k,
        )?;
        let contributions = values
            .into_iter()
            .map(|(key, value)| {
                WindowBoundaryContribution::try_new(
                    &spec,
                    &id(&key),
                    value,
                    evidence(
                        &interval,
                        WindowEvidenceRole::IntegratedTransfer {
                            contribution_id: id(&key),
                            local_port_id: manifest
                                .iter()
                                .find(|m| m.contribution_id() == &id(&key))
                                .expect("declared row")
                                .local_port_id()
                                .clone(),
                            energy_reference: reference.clone(),
                        },
                        &format!("integration/{key}"),
                    ),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(ModuleWindowAudit::Available(Box::new(
            WindowAuditReport::audit(spec, initial, final_snapshot, contributions)?,
        )))
    }
}

fn evidence(
    interval: &AccountingWindowInterval,
    role: WindowEvidenceRole,
    key: &str,
) -> WindowEvidenceRef {
    WindowEvidenceRef::new(
        id(&format!("{}/{key}", interval.id().as_str())),
        id("fs-couple/accepted-step-adapter"),
        id(&format!("source-account/{key}")),
        role,
        interval.clone(),
    )
}
fn boundary_entry(
    port: &PortSchema,
    key: &str,
    interval: &AccountingWindowInterval,
) -> Result<WindowManifestEntry, CoupleError> {
    WindowManifestEntry::external_power_port(
        id(key),
        id(&format!("exchange/{key}")),
        id(&format!("environment/{key}")),
        port,
        AccountingBoundary::new(
            id(&format!("boundary/{key}")),
            port.coordinates().clone(),
            BoundaryTreatment::ExternalReservoirExchange,
        ),
        evidence(
            interval,
            WindowEvidenceRole::ManifestBinding {
                contribution_id: id(key),
                local_port_id: port.id().clone(),
            },
            &format!("binding/{key}"),
        ),
    )
}
fn transfer(
    energy: f64,
    temperature: Option<f64>,
    interval: &AccountingWindowInterval,
    key: &str,
) -> Result<IntegratedWindowTransfer, CoupleError> {
    let (heat, temperature_reference) = match temperature {
        Some(t) => {
            if !t.is_finite() || t <= 0.0 {
                return Err(CoupleError::NonPositiveBoundaryTemperature);
            }
            (
                energy / t,
                BoundaryTemperatureReference::profile(
                    id(&format!("midpoint-temperature/{:016x}", t.to_bits())),
                    evidence(
                        interval,
                        WindowEvidenceRole::BoundaryTemperature {
                            contribution_id: id(key),
                        },
                        &format!("temperature/{key}"),
                    ),
                ),
            )
        }
        None => (
            0.0,
            BoundaryTemperatureReference::not_applicable(id("mechanical/no-heat-entropy")),
        ),
    };
    IntegratedWindowTransfer::try_new(
        Energy::new(energy),
        Mass::new(0.0),
        [],
        ElectricCharge::new(0.0),
        BoundaryEntropyBreakdown::try_new(
            Entropy::new(0.0),
            Entropy::new(0.0),
            Entropy::new(heat),
            Entropy::new(0.0),
            temperature_reference,
        )?,
        id("accepted-step/midpoint-heat-entropy"),
    )
}

macro_rules! impl_ported_module {
    ($owner:ty, $field:ident) => {
        impl fs_ports::PortedModule for $owner {
            fn port_schemas(&self) -> &[fs_ports::PortSchema] {
                &self.$field.schemas
            }
            fn port_coordinate_charts(&self) -> &[fs_ports::PortCoordinateChart] {
                &self.$field.charts
            }
            fn junction_schemas(&self) -> &[fs_ports::JunctionSchema] {
                &self.$field.junctions
            }
            fn last_window_duration_s(&self) -> Option<f64> {
                self.$field.duration_s()
            }
            fn last_window_transfers(&self) -> Option<&[fs_ports::IntegratedWindowTransfer]> {
                self.$field.transfers()
            }
            fn last_window_transfer_error(&self) -> Option<&fs_ports::CoupleError> {
                self.$field.transfer_error()
            }
            fn audit_last_window(
                &self,
            ) -> Result<fs_ports::ModuleWindowAudit, fs_ports::CoupleError> {
                self.$field.audit()
            }
        }
    };
}
pub(crate) use impl_ported_module;

pub(crate) fn junction(
    owner: &str,
    name: &str,
    count: usize,
    kind: PortKind,
    common_flow: bool,
) -> JunctionSchema {
    let ports = (0..count)
        .map(|i| schema(owner, &format!("junction/{name}/{i}"), kind))
        .collect();
    let law = if common_flow {
        JunctionLaw::CommonFlow {
            effort_signs: vec![1; count],
        }
    } else {
        JunctionLaw::CommonEffort {
            flow_signs: vec![1; count],
        }
    };
    JunctionSchema::try_new(id(&format!("{owner}/junction/{name}")), ports, law)
        .expect("admitted owner count, canonical conjugate charts and unit signs")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_observation_refusal_is_visible_and_next_acceptance_replaces_it() {
        let mut ports = ModulePorts::mechanical("observation-test", 1);
        let account = StepAccount {
            tick: 1,
            duration_s: 0.01,
            storage_change_j: Some(0.0),
            loss_j: Some(0.0),
            included_work_j: 0.0,
            entropy_change_j_k: 0.0,
            tolerance_j: Some(1e-12),
            unavailable: None,
        };
        let error = CoupleError::ModuleObservationFailed {
            owner: "observation-test".into(),
            reason: "independent loss integral refused".into(),
        };
        ports.accept_result(account.clone(), Err(error.clone()));
        assert_eq!(ports.duration_s(), Some(0.01));
        assert!(ports.transfers().is_none());
        assert_eq!(ports.transfer_error(), Some(&error));
        assert_eq!(ports.audit(), Err(error));
        ports.accept_result(
            StepAccount { tick: 2, ..account },
            Ok(vec![PortTransfer::work(0.0, true)]),
        );
        assert_eq!(ports.next_tick(), 3);
        assert!(ports.transfer_error().is_none());
        assert_eq!(ports.transfers().unwrap()[0].energy_into().value(), 0.0);
        let ModuleWindowAudit::Available(report) = ports.audit().unwrap() else {
            panic!("replacement accepted account missing")
        };
        assert!(report.violations().is_empty());
    }
}
