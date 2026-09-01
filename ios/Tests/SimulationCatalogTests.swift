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
            XCTAssertEqual(error.localizedDescription, "Invalid PCM payload: samples must be finite normalized PCM")
        }
    }

    func testStartupPrefillCallbackDefersPlaybackUntilFourBuffersAreQueued() {
        var schedule = PCMStreamSchedule(plan: .reedDemo)
        for _ in 0..<(PCMStreamPlan.reedDemo.maximumQueuedBlocks - 1) {
            XCTAssertTrue(schedule.reserveNextBlock())
            XCTAssertFalse(schedule.isReadyToStartPlayback)
        }
        XCTAssertEqual(schedule.scheduledBlocks, 3)
        XCTAssertEqual(schedule.queuedBlocks, 3)

        XCTAssertTrue(schedule.reserveNextBlock())
        XCTAssertEqual(schedule.queuedBlocks, 4)
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
        XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
        XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
        XCTAssertEqual(schedule.didDrainBlock(), .awaitingRefill)
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
