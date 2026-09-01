import XCTest
@testable import FrankenSim
import FrankenSimCore

final class SimulationCatalogTests: XCTestCase {
    func testCompleteWebsiteAndCampaignInventoryIsUnique() {
        XCTAssertEqual(SimulationCatalog.all.count, 43)
        XCTAssertEqual(Set(SimulationCatalog.all.map(\.id)).count, 43)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .foundations }.count, 10)
        XCTAssertEqual(SimulationCatalog.all.filter { $0.tier == .frontier }.count, 10)
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
}
