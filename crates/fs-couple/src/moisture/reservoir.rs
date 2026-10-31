//! Small field adapter; finite vapor and solution inventories keep their owners.

use super::{
    atmosphere::DryBasisSorption,
    binding::MoistureField,
    cartesian::{MoistureCartesian, MoistureCartesianReport},
    surface::ResolvedMoistureFace,
    *,
};

/// Compensated signed water inventories, including opposing specimen transfers.
/// A small net reservoir change must not erase represented local exchanges.
pub(super) fn water_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

/// Only the two original physical fields implement this private adapter.
pub(super) trait ReservoirField<const FACES: usize>: Sized {
    type Report;
    type Controls: Copy;
    fn inventory(&self) -> MoistureField<'_>;
    fn advance_resolved(
        &self,
        dt_s: f64,
        boundary: [ResolvedMoistureFace; FACES],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError>;

    fn advance_activity(
        &self,
        dt_s: f64,
        laws: &[Option<DryBasisSorption>; FACES],
        temperature: f64,
        activity: f64,
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        let mut boundary = std::array::from_fn(|_| ResolvedMoistureFace::default());
        for (side, law) in laws.iter().enumerate() {
            poll(cx)?;
            if let Some(law) = law {
                boundary[side] = law.resolve(temperature, activity)?;
            }
        }
        self.advance_resolved(dt_s, boundary, controls, cx)
    }

    /// Use represented local changes, not subtraction of large stored totals.
    fn water_change_to(&self, next: &Self) -> f64 {
        let old = self.inventory();
        water_sum(
            next.inventory()
                .moisture
                .iter()
                .zip(old.moisture)
                .map(|(new, before)| (new - before) * old.dry_mass_per_cell_kg),
        )
    }
}

impl ReservoirField<2> for MoistureSlab {
    type Report = MoistureStepReport;
    type Controls = MoistureStepConfig;
    fn inventory(&self) -> MoistureField<'_> {
        self.field()
    }
    fn advance_resolved(
        &self,
        dt_s: f64,
        boundary: [ResolvedMoistureFace; 2],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        MoistureSlab::advance_resolved(self, dt_s, boundary, controls, cx)
    }
}

impl ReservoirField<6> for MoistureCartesian {
    type Report = MoistureCartesianReport;
    type Controls = (MoistureStepConfig, MoistureNonlinearConfig);
    fn inventory(&self) -> MoistureField<'_> {
        self.field()
    }
    fn advance_resolved(
        &self,
        dt_s: f64,
        boundary: [ResolvedMoistureFace; 6],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        MoistureCartesian::advance_resolved(self, dt_s, boundary, controls.0, controls.1, cx)
    }
}

pub(super) struct ReservoirParticipant<'a, F: ReservoirField<N>, const N: usize> {
    pub id: u64,
    pub field: &'a F,
    pub faces: &'a [Option<DryBasisSorption>; N],
    pub transport: F::Controls,
}
