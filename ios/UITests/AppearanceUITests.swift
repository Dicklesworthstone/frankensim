import XCTest

final class FrankenSimAppearanceUITests: XCTestCase {
    override func setUpWithError() throws {
        continueAfterFailure = false
    }

    func testAppearanceTogglePersistsLightModeAcrossLaunches() throws {
        let app = XCUIApplication()
        app.launch()

        let toggle = app.buttons["appearance-toggle"]
        XCTAssertTrue(toggle.waitForExistence(timeout: 12))
        XCTAssertTrue(
            ["Switch to light mode", "Switch to dark mode"].contains(toggle.label),
            "Appearance control exposed an unexpected state: \(toggle.label)"
        )

        if toggle.label == "Switch to dark mode" {
            toggle.tap()
            XCTAssertEqual(toggle.label, "Switch to light mode")
        }

        toggle.tap()
        XCTAssertEqual(toggle.label, "Switch to dark mode")
        keepScreenshot(of: app, named: "Remembered light appearance")

        app.terminate()
        app.launch()

        let relaunchedToggle = app.buttons["appearance-toggle"]
        XCTAssertTrue(relaunchedToggle.waitForExistence(timeout: 12))
        XCTAssertEqual(relaunchedToggle.label, "Switch to dark mode")
    }

    func testPhoneStudioExposesAnExplicitRunAndTheFullCatalog() throws {
        let app = XCUIApplication()
        app.launchEnvironment["FSIM_INITIAL_EXPERIMENT"] = "13"
        app.launch()

        let chooser = app.buttons["compact-kernel-chooser"]
        XCTAssertTrue(chooser.waitForExistence(timeout: 12))
        XCTAssertTrue(app.otherElements["Simulation canvas ready"].exists)
        XCTAssertFalse(app.otherElements["Rendered native simulation result"].exists)

        let run = app.buttons["compact-run-button"]
        XCTAssertTrue(run.waitForExistence(timeout: 5))
        XCTAssertTrue(run.isHittable)
        XCTAssertEqual(run.label, "Run Lorenz attractor")
        XCTAssertEqual(run.value as? String, "Ready")
        run.tap()
        XCTAssertTrue(app.otherElements["Rendered native simulation result"].waitForExistence(timeout: 12))
        XCTAssertEqual(run.label, "Run Lorenz attractor")
        XCTAssertEqual(run.value as? String, "Completed")
        XCTAssertTrue(run.isHittable)
        keepScreenshot(of: app, named: "Phone studio after explicit run")

        chooser.tap()
        XCTAssertTrue(app.navigationBars["Choose a simulation"].waitForExistence(timeout: 5))
        XCTAssertTrue(app.staticTexts["Heat diffusion"].exists)
        XCTAssertTrue(app.searchFields.firstMatch.exists)
        keepScreenshot(of: app, named: "Phone simulation catalog")
    }

    private func keepScreenshot(of app: XCUIApplication, named name: String) {
        let attachment = XCTAttachment(screenshot: app.screenshot())
        attachment.name = name
        attachment.lifetime = .keepAlways
        add(attachment)
    }
}
