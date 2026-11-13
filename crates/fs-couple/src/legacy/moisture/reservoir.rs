//! Small field adapter; finite vapor and solution inventories keep their owners.

use super::{
    atmosphere::DryBasisSorption,
    binding::MoistureField,
    cartesian::{MoistureCartesian, MoistureCartesianReport},
    radial::MoistureRadial,
    sdirk::StageLoad,
    surface::ResolvedMoistureFace,
    *,
};
use std::borrow::Cow;

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

/// Native fields and borrowed implicit stages share this private adapter.
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
                .enumerate()
                .map(|(i, (new, before))| (new - before) * old.dry_mass_at(i)),
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

impl ReservoirField<2> for MoistureRadial {
    type Report = MoistureStepReport;
    type Controls = (MoistureStepConfig, MoistureNonlinearConfig);

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
        self.advance_stage(dt_s, dt_s, boundary, controls.0, controls.1, None, cx)
    }
}

pub(super) struct ReservoirParticipant<'a, F: ReservoirField<N>, const N: usize> {
    pub id: u64,
    pub field: &'a F,
    pub faces: &'a [Option<DryBasisSorption>; N],
    pub transport: F::Controls,
}

// Borrow the original storage/load during scalar reservoir trials. Successful
// proposals own the native field; source fields and loads are not cloned per trial.
pub(super) struct SlabSdirkStage<'a> {
    field: Cow<'a, MoistureSlab>,
    known: &'a StageLoad,
    implicit_dt_s: f64,
}

impl<'a> SlabSdirkStage<'a> {
    pub(super) fn new(
        field: &'a MoistureSlab,
        known: &'a StageLoad,
        implicit_dt_s: f64,
    ) -> Self {
        Self { field: Cow::Borrowed(field), known, implicit_dt_s }
    }

    pub(super) fn into_field(self) -> MoistureSlab {
        self.field.into_owned()
    }
}

impl ReservoirField<2> for SlabSdirkStage<'_> {
    type Report = MoistureStepReport;
    type Controls = MoistureStepConfig;

    fn inventory(&self) -> MoistureField<'_> {
        self.field.field()
    }

    fn advance_resolved(
        &self,
        dt_s: f64,
        boundary: [ResolvedMoistureFace; 2],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        let (field, report) = self.field.advance_stage(
            dt_s, self.implicit_dt_s, boundary, controls, Some(self.known), cx,
        )?;
        Ok((Self { field: Cow::Owned(field), known: self.known,
            implicit_dt_s: self.implicit_dt_s }, report))
    }
}

pub(super) struct CartesianSdirkStage<'a> {
    field: Cow<'a, MoistureCartesian>,
    known: &'a StageLoad<6>,
    implicit_dt_s: f64,
}

impl<'a> CartesianSdirkStage<'a> {
    pub(super) fn new(
        field: &'a MoistureCartesian,
        known: &'a StageLoad<6>,
        implicit_dt_s: f64,
    ) -> Self {
        Self {
            field: Cow::Borrowed(field),
            known,
            implicit_dt_s,
        }
    }

    pub(super) fn into_field(self) -> MoistureCartesian {
        self.field.into_owned()
    }
}

impl ReservoirField<6> for CartesianSdirkStage<'_> {
    type Report = MoistureCartesianReport;
    type Controls = (MoistureStepConfig, MoistureNonlinearConfig);

    fn inventory(&self) -> MoistureField<'_> {
        self.field.field()
    }

    fn advance_resolved(
        &self,
        dt_s: f64,
        boundaries: [ResolvedMoistureFace; 6],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        let (field, report) = self.field.advance_stage(
            dt_s,
            self.implicit_dt_s,
            boundaries,
            controls.0,
            controls.1,
            Some(self.known),
            cx,
        )?;
        Ok((
            Self {
                field: Cow::Owned(field),
                known: self.known,
                implicit_dt_s: self.implicit_dt_s,
            },
            report,
        ))
    }
}

/// Borrow original curved-shell storage during finite-reservoir stage trials.
pub(super) struct RadialSdirkStage<'a> {
    field: Cow<'a, MoistureRadial>,
    known: &'a StageLoad,
    implicit_dt_s: f64,
}

impl<'a> RadialSdirkStage<'a> {
    pub(super) fn new(field: &'a MoistureRadial, known: &'a StageLoad, implicit_dt_s: f64) -> Self {
        Self {
            field: Cow::Borrowed(field),
            known,
            implicit_dt_s,
        }
    }

    pub(super) fn into_field(self) -> MoistureRadial {
        self.field.into_owned()
    }
}

impl ReservoirField<2> for RadialSdirkStage<'_> {
    type Report = MoistureStepReport;
    type Controls = (MoistureStepConfig, MoistureNonlinearConfig);

    fn inventory(&self) -> MoistureField<'_> {
        self.field.field()
    }

    fn advance_resolved(
        &self,
        dt_s: f64,
        boundary: [ResolvedMoistureFace; 2],
        controls: Self::Controls,
        cx: &Cx<'_>,
    ) -> Result<(Self, Self::Report), MoistureError> {
        let (field, report) = self.field.advance_stage(
            dt_s,
            self.implicit_dt_s,
            boundary,
            controls.0,
            controls.1,
            Some(self.known),
            cx,
        )?;
        Ok((
            Self {
                field: Cow::Owned(field),
                known: self.known,
                implicit_dt_s: self.implicit_dt_s,
            },
            report,
        ))
    }
}
