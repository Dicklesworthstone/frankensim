//! One molecular dispatch for complete network compositions.

use super::super::MixtureNetwork1d;
use super::super::dissipation::{
    AdaptiveNetworkDissipationControl, AdaptiveNetworkDissipationReport,
    NetworkDissipationOperator, NetworkDissipationReport, SplitNetworkDissipationReport,
    StabilizedNetworkDissipationReport,
};
use crate::multicomponent::Cell;
use crate::multicomponent::dissipation::{
    DissipationControl, DissipationError, StabilizedDissipationControl, TransportLaw,
};
use fs_qty::Time;

/// Accepted RKL2 intervals on both sides of the original network composition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitStabilizedNetworkDissipationReport {
    pub first: StabilizedNetworkDissipationReport,
    pub second: StabilizedNetworkDissipationReport,
}

/// Accepted diffusion error estimates on both sides of the network composition.
/// These estimate molecular time error, not the splitting, chemistry or wave error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SplitAdaptiveNetworkDissipationReport {
    pub first: AdaptiveNetworkDissipationReport,
    pub second: AdaptiveNetworkDissipationReport,
}

#[derive(Default)]
pub(in crate::multicomponent::network) struct MolecularDiagnostics {
    pub stabilized: Option<SplitStabilizedNetworkDissipationReport>,
    pub adaptive: Option<SplitAdaptiveNetworkDissipationReport>,
}

pub(in crate::multicomponent::network) enum MolecularMethod<'a> {
    Explicit(DissipationControl),
    Stabilized {
        controls: StabilizedDissipationControl,
        /// Borrowed by the complete composition, never reset by a half or retry.
        evaluations: &'a mut usize,
    },
    Adaptive {
        controls: AdaptiveNetworkDissipationControl,
        /// Includes all temporal probes, both halves and discarded compositions.
        evaluations: &'a mut usize,
    },
}

pub(in crate::multicomponent::network) struct MolecularHalf {
    physical: NetworkDissipationReport,
    stabilized: Option<StabilizedNetworkDissipationReport>,
    adaptive: Option<AdaptiveNetworkDissipationReport>,
}

impl MolecularMethod<'_> {
    #[allow(clippy::too_many_arguments)]
    pub(in crate::multicomponent::network) fn advance<L: TransportLaw>(
        &mut self,
        operator: &mut NetworkDissipationOperator<L>,
        network: &mut MixtureNetwork1d,
        duration: Time,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        moving_volumes: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularHalf, DissipationError> {
        match self {
            Self::Explicit(controls) => {
                let physical = if moving_volumes {
                    operator.advance_for_moving_trial(
                        network, duration, *controls, reservoirs, cancelled,
                    )?
                } else {
                    operator
                        .advance_for_trial(network, duration, *controls, reservoirs, cancelled)?
                };
                Ok(MolecularHalf {
                    physical,
                    stabilized: None,
                    adaptive: None,
                })
            }
            Self::Stabilized {
                controls,
                evaluations,
            } => {
                let report = if moving_volumes {
                    operator.advance_stabilized_for_moving_trial(
                        network,
                        duration,
                        *controls,
                        reservoirs,
                        evaluations,
                        cancelled,
                    )?
                } else {
                    operator.advance_stabilized_for_trial(
                        network,
                        duration,
                        *controls,
                        reservoirs,
                        evaluations,
                        cancelled,
                    )?
                };
                Ok(MolecularHalf {
                    physical: report.transport,
                    stabilized: Some(report),
                    adaptive: None,
                })
            }
            Self::Adaptive {
                controls,
                evaluations,
            } => {
                let report = if moving_volumes {
                    operator.advance_adaptive_for_moving_trial(
                        network,
                        duration,
                        *controls,
                        reservoirs,
                        evaluations,
                        cancelled,
                    )?
                } else {
                    operator.advance_adaptive_for_trial(
                        network,
                        duration,
                        *controls,
                        reservoirs,
                        evaluations,
                        cancelled,
                    )?
                };
                Ok(MolecularHalf {
                    physical: report.stabilized.transport,
                    stabilized: Some(report.stabilized),
                    adaptive: Some(report),
                })
            }
        }
    }
}

pub(in crate::multicomponent::network) fn split_halves(
    first: MolecularHalf,
    second: MolecularHalf,
) -> (SplitNetworkDissipationReport, MolecularDiagnostics) {
    (
        SplitNetworkDissipationReport {
            first: first.physical,
            second: second.physical,
        },
        MolecularDiagnostics {
            stabilized: first
                .stabilized
                .zip(second.stabilized)
                .map(|(first, second)| SplitStabilizedNetworkDissipationReport { first, second }),
            adaptive: first
                .adaptive
                .zip(second.adaptive)
                .map(|(first, second)| SplitAdaptiveNetworkDissipationReport { first, second }),
        },
    )
}
