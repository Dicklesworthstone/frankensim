import Foundation
import XCTest
@testable import FrankenSim
import FrankenSimCore

final class SimulationCatalogTests: XCTestCase {
    func testCompleteWebsiteAndCampaignInventoryIsUnique() {
        XCTAssertEqual(SimulationCatalog.all.count, 44)
        XCTAssertEqual(Set(SimulationCatalog.all.map(\.id)).count, 44)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .foundations }.count, 10)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .frontier }.count, 11)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .deep }.count, 10)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .campaigns }.count, 10)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .flagships }.count, 3)
    }

    func testEveryExperimentStatesANoClaimBoundary() {
        XCTAssertTrue(SimulationCatalog.all.allSatisfy { !$0.noClaim.isEmpty && !$0.kernel.isEmpty })
    }

    func testCatalogFiltersExposeEveryTierWithoutDroppingEntries() {
        let filteredInventory = ExperimentTier.allCases.flatMap { tier in
            SimulationCatalog.filtered(query: "", tier: tier)
        }

        XCTAssertEqual(filteredInventory.count, SimulationCatalog.all.count)
        XCTAssertEqual(Set(filteredInventory.map(\.id)), Set(SimulationCatalog.all.map(\.id)))
        XCTAssertEqual(SimulationCatalog.filtered(query: "", tier: .flagships).map(\.id), [40, 41, 42])
    }

    func testCatalogSearchCoversMethodsEvidenceAndExplanations() {
        XCTAssertEqual(SimulationCatalog.filtered(query: "fs-couple").map(\.id), [43])
        XCTAssertEqual(SimulationCatalog.filtered(query: "Verified bound").map(\.id), [3, 9, 17, 29])
        XCTAssertTrue(SimulationCatalog.filtered(query: "buoyant smoke").contains { $0.id == 15 })
        let tierSearch = SimulationCatalog.filtered(query: "Tier III")
        XCTAssertEqual(tierSearch.count, 10)
        XCTAssertTrue(tierSearch.allSatisfy { $0.tier == .deep })
        XCTAssertTrue(SimulationCatalog.filtered(query: "does not exist").isEmpty)
    }

    @MainActor
    func testSelectingAnExperimentWaitsForAnExplicitRun() throws {
        let model = SimulationStudioModel()
        let experiment = try XCTUnwrap(SimulationCatalog.all.first { $0.id == 17 })

        model.select(experiment)

        XCTAssertEqual(model.selection.id, experiment.id)
        XCTAssertNil(model.result)
        XCTAssertFalse(model.isRunning)
    }

    func testNativeImagePacketHasTheSixValueABIHeader() {
        let header = runAndCopyHeader(id: 17)

        XCTAssertEqual(header, [1, 17, 1, 96, 96, 1])
    }

    func testNativeSeriesPacketHasTheSixValueABIHeader() {
        let header = runAndCopyHeader(id: 1)

        XCTAssertEqual(header, [1, 1, 0, 96, 1, 1])
    }

    func testNativeScalarPacketHasTheSixValueABIHeader() {
        let header = runAndCopyHeader(id: 3)

        XCTAssertEqual(header, [1, 3, 0, 14, 1, 1])
    }

    func testNativeReedPacketIsTypedMonoPCM() {
        let header = runAndCopyHeader(id: 43)

        XCTAssertEqual(header, [1, 43, 6, 480, 1, 48_000])
    }

    func testPCMBlockAcceptsOnlyTheDeclaredMonoFormatAndFiniteSamples() throws {
        let block = try PCMBlock(result: pcmResult(values: [0, -0.5, 1]))
        XCTAssertEqual(block.samples, [0, -0.5, 1])

        XCTAssertThrowsError(try PCMBlock(result: pcmResult(channels: 2))) { error in
            XCTAssertEqual(error.localizedDescription, "Unsupported PCM format: 2 channel(s) at 48000 Hz.")
        }
        XCTAssertThrowsError(try PCMBlock(result: pcmResult(values: [Double.nan]))) { error in
            XCTAssertEqual(error.localizedDescription, "Invalid PCM payload: samples must be finite normalized PCM.")
        }
    }

    func testStartupPrefillCallbackDefersPlaybackUntilBoundedLeadIsQueued() {
        var schedule = PCMStreamSchedule(plan: .reedDemo)
        for _ in 0..<(PCMStreamPlan.reedDemo.maximumQueuedBlocks - 1) {
            XCTAssertTrue(schedule.reserveNextBlock())
            XCTAssertFalse(schedule.isReadyToStartPlayback)
        }
        XCTAssertEqual(
            schedule.scheduledBlocks,
            PCMStreamPlan.reedDemo.maximumQueuedBlocks - 1
        )
        XCTAssertEqual(
            schedule.queuedBlocks,
            PCMStreamPlan.reedDemo.maximumQueuedBlocks - 1
        )

        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertEqual(schedule.queuedBlocks, PCMStreamPlan.reedDemo.maximumQueuedBlocks)
        XCTAssertTrue(schedule.isReadyToStartPlayback)
        schedule.markPlaybackStarted()
        XCTAssertTrue(schedule.playbackStarted)
        XCTAssertFalse(schedule.isReadyToStartPlayback)
    }

    func testRefillAndDrainCallbacksCompleteFiniteStreamWithoutADevice() {
        var schedule = PCMStreamSchedule(
            plan: PCMStreamPlan(totalBlocks: 3, maximumQueuedBlocks: 2)
        )
        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertTrue(schedule.isReadyToStartPlayback)
        schedule.markPlaybackStarted()
        XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
        XCTAssertEqual(schedule.didDrainBlock(), .complete)
        XCTAssertEqual(schedule.scheduledBlocks, 3)
        XCTAssertEqual(schedule.drainedBlocks, 3)
    }

    func testStopCallbackDrainsAllProducedBuffersBeforeStopping() {
        var schedule = PCMStreamSchedule(plan: .reedDemo)
        for _ in 0..<PCMStreamPlan.reedDemo.maximumQueuedBlocks {
            XCTAssertTrue(schedule.reserveNextBlock())
        }
        schedule.markPlaybackStarted()
        schedule.requestStopAfterDrain()
        for _ in 1..<PCMStreamPlan.reedDemo.maximumQueuedBlocks {
            XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
        }
        XCTAssertEqual(schedule.didDrainBlock(), .stopped)
        XCTAssertFalse(schedule.reserveNextBlock())
    }

    func testPreplayStopTerminatesPartialPrefillWithoutADevice() {
        var schedule = PCMStreamSchedule(plan: .reedDemo)
        for _ in 0..<(PCMStreamPlan.reedDemo.maximumQueuedBlocks - 1) {
            XCTAssertTrue(schedule.reserveNextBlock())
        }

        XCTAssertTrue(schedule.stopsImmediately)
        schedule.requestStopAfterDrain()
        XCTAssertTrue(schedule.stopping)
        XCTAssertFalse(schedule.reserveNextBlock())
    }

    func testEmptyQueueBeforeFiniteProductionCompletesIsStarvation() {
        var schedule = PCMStreamSchedule(
            plan: PCMStreamPlan(totalBlocks: 3, maximumQueuedBlocks: 1)
        )
        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertEqual(schedule.didDrainBlock(), .starved)
    }

    func testStaleCompletionCallbackIsRejectedAfterStreamRestart() {
        var generation = PCMStreamGeneration()
        let staleStreamID = generation.begin()
        let replacementStreamID = generation.begin()

        XCTAssertFalse(generation.accepts(staleStreamID))
        XCTAssertTrue(generation.accepts(replacementStreamID))
    }

    func testNativeBridgeRefusesAnUnknownCatalogID() async {
        let unknown = SimulationExperiment(
            id: .max,
            name: "unknown",
            subtitle: "",
            explanation: "",
            tier: .foundations,
            symbol: "",
            accent: .cyan,
            evidence: .estimated,
            noClaim: "",
            kernel: ""
        )

        do {
            _ = try await NativeSimulationKernel().run(experiment: unknown, quality: 0, seed: 1)
            XCTFail("the bridge must refuse an unknown catalog id")
        } catch let error as KernelError {
            XCTAssertEqual(error.errorDescription, Optional("The native kernel refused this run (code 1)."))
            XCTAssertEqual(frankensim_apple_result_len(), 0)
        } catch {
            XCTFail("unexpected error: \(error)")
        }
    }

    func testRunSnapshotExportsCompletePayloadAndAuthorityBoundary() throws {
        let experiment = try XCTUnwrap(SimulationCatalog.all.first { $0.id == 13 })
        let values = [1.25, .nan, .infinity, -.infinity]
        let result = SimulationResult(
            experimentID: experiment.id,
            shape: .signal,
            width: values.count,
            height: 1,
            frames: 1,
            values: values,
            finiteMinimum: 1.25,
            finiteMaximum: 1.25,
            quality: 0.55,
            seed: 0x5EED,
            elapsed: .milliseconds(125)
        )

        let snapshot = try SimulationRunSnapshot(
            result: result,
            experiment: experiment,
            createdAt: Date(timeIntervalSince1970: 2)
        )
        let json = try XCTUnwrap(
            JSONSerialization.jsonObject(with: snapshot.encoded()) as? [String: Any]
        )
        let payload = try XCTUnwrap(json["payload"] as? [Any])
        let run = try XCTUnwrap(json["run"] as? [String: Any])

        XCTAssertEqual(json["schema"] as? String, SimulationRunSnapshot.schema)
        XCTAssertEqual(json["authorityBoundary"] as? String, SimulationRunSnapshot.boundary)
        XCTAssertEqual(json["noClaim"] as? String, experiment.noClaim)
        XCTAssertEqual(payload.count, values.count)
        XCTAssertEqual(payload[0] as? Double, 1.25)
        XCTAssertEqual(payload[1] as? String, "nan")
        XCTAssertEqual(payload[2] as? String, "+infinity")
        XCTAssertEqual(payload[3] as? String, "-infinity")
        XCTAssertEqual(run["qualityTier"] as? String, "balanced")
        XCTAssertEqual(run["elapsedMilliseconds"] as? Double, 125)
        XCTAssertTrue(snapshot.filename.hasSuffix("-5eed.json"))
    }

    func testRunSnapshotRefusesMismatchedExperimentIdentity() throws {
        let experiment = try XCTUnwrap(SimulationCatalog.all.first { $0.id == 13 })
        let result = SimulationResult(
            experimentID: 12,
            shape: .signal,
            width: 1,
            height: 1,
            frames: 1,
            values: [1],
            finiteMinimum: 1,
            finiteMaximum: 1,
            quality: 0.12,
            seed: 1,
            elapsed: .zero
        )

        XCTAssertThrowsError(
            try SimulationRunSnapshot(result: result, experiment: experiment)
        ) { error in
            XCTAssertEqual(
                error.localizedDescription,
                "The result does not belong to the selected experiment."
            )
        }
    }

    func testResultPlaybackAdvancesLoopsAndPausesOnAnExactFrame() {
        let start = Date(timeIntervalSinceReferenceDate: 100)
        var playback = ResultPlaybackState(frameCount: 5, date: start)

        XCTAssertTrue(playback.isPlaying)
        XCTAssertEqual(playback.displayedFrame(at: start), 0)
        XCTAssertEqual(playback.displayedFrame(at: start.addingTimeInterval(0.25)), 2)
        XCTAssertEqual(playback.displayedFrame(at: start.addingTimeInterval(1)), 4)

        playback.pause(at: start.addingTimeInterval(0.25))
        XCTAssertFalse(playback.isPlaying)
        XCTAssertEqual(playback.displayedFrame(at: start.addingTimeInterval(30)), 2)
    }

    func testResultPlaybackSeekRestartRateAndReduceMotionTransitions() {
        let start = Date(timeIntervalSinceReferenceDate: 200)
        var playback = ResultPlaybackState(frameCount: 6, reduceMotion: true, date: start)

        XCTAssertFalse(playback.isPlaying)
        XCTAssertEqual(playback.displayedFrame(at: start), 5)

        playback.seek(to: 99, at: start)
        XCTAssertEqual(playback.displayedFrame(at: start), 5)
        playback.seek(to: 2, at: start)
        XCTAssertEqual(playback.displayedFrame(at: start), 2)

        playback.toggle(at: start, reduceMotion: true)
        XCTAssertFalse(playback.isPlaying)
        playback.restart(at: start, reduceMotion: false)
        XCTAssertTrue(playback.isPlaying)
        XCTAssertEqual(playback.displayedFrame(at: start), 0)

        playback.setRate(2, at: start)
        XCTAssertEqual(playback.rate, 2)
        XCTAssertEqual(playback.displayedFrame(at: start.addingTimeInterval(0.25)), 4)
        playback.setRate(3, at: start.addingTimeInterval(0.25))
        XCTAssertEqual(playback.rate, 2, "unsupported presentation rates must not be admitted")

        playback.honorReduceMotion(true, at: start.addingTimeInterval(0.25))
        XCTAssertFalse(playback.isPlaying)
        XCTAssertEqual(playback.displayedFrame(at: start.addingTimeInterval(10)), 4)
    }

    private func runAndCopyHeader(id: UInt32) -> [Double] {
        let reported = frankensim_apple_run(id, 0, 0x5EED)
        XCTAssertEqual(frankensim_apple_last_error(), 0)
        XCTAssertEqual(reported, frankensim_apple_result_len())
        XCTAssertGreaterThanOrEqual(reported, 6)

        let header = (0..<6).map { frankensim_apple_result_value(UInt64($0)) }
        XCTAssertTrue(header.allSatisfy(\.isFinite))
        return header
    }

    private func pcmResult(
        values: [Double] = [0],
        channels: Int = 1,
        sampleRate: Int = 48_000
    ) -> SimulationResult {
        SimulationResult(
            experimentID: 43,
            shape: .pcm,
            width: values.count,
            height: channels,
            frames: sampleRate,
            values: values,
            finiteMinimum: values.filter(\.isFinite).min() ?? 0,
            finiteMaximum: values.filter(\.isFinite).max() ?? 0,
            quality: 0,
            seed: 0,
            elapsed: .zero
        )
    }
}
