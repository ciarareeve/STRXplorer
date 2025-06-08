#!/usr/bin/env python3
"""
STRXplorer Endpoint Testing Script
Automatically tests all endpoints and logs results
Updated to use correct database data
"""
import requests
import json
import time
import sys
from datetime import datetime
from typing import Dict, List, Tuple
import sqlite3
import os


class EndpointTester:
    def __init__(self, base_url: str = "http://localhost:5000"):
        self.base_url = base_url.rstrip("/")
        self.results = []
        self.errors = []
        self.start_time = datetime.now()

    def log(self, message: str, level: str = "INFO"):
        """Log a message with timestamp"""
        timestamp = datetime.now().strftime("%H:%M:%S")
        print(f"[{timestamp}] {level}: {message}")

    def test_endpoint(
        self,
        endpoint: str,
        method: str = "GET",
        expected_status: int = 200,
        params: Dict = None,
        description: str = None,
    ) -> bool:
        """Test a single endpoint"""
        url = f"{self.base_url}{endpoint}"
        test_name = description or f"{method} {endpoint}"

        try:
            if method == "GET":
                response = requests.get(url, params=params, timeout=10)
            elif method == "POST":
                response = requests.post(url, data=params, timeout=10)
            else:
                raise ValueError(f"Unsupported method: {method}")

            success = response.status_code == expected_status

            result = {
                "test": test_name,
                "endpoint": endpoint,
                "method": method,
                "status_code": response.status_code,
                "expected_status": expected_status,
                "success": success,
                "response_time_ms": response.elapsed.total_seconds() * 1000,
                "content_length": len(response.content),
                "content_type": response.headers.get("content-type", "unknown"),
            }

            if success:
                self.log(
                    f"{test_name} - {response.status_code} ({result['response_time_ms']:.0f}ms)"
                )
            else:
                self.log(
                    f"❌ {test_name} - Expected {expected_status}, got {response.status_code}",
                    "ERROR",
                )
                self.errors.append(
                    {
                        "test": test_name,
                        "error": f"Status code mismatch: {response.status_code} != {expected_status}",
                        "response_text": response.text[:500]
                        if response.text
                        else "No response text",
                    }
                )

            self.results.append(result)
            return success

        except requests.exceptions.RequestException as e:
            self.log(f"❌ {test_name} - Connection error: {str(e)}", "ERROR")
            self.errors.append(
                {
                    "test": test_name,
                    "error": f"Connection error: {str(e)}",
                    "response_text": None,
                }
            )
            return False
        except Exception as e:
            self.log(f"❌ {test_name} - Unexpected error: {str(e)}", "ERROR")
            self.errors.append(
                {
                    "test": test_name,
                    "error": f"Unexpected error: {str(e)}",
                    "response_text": None,
                }
            )
            return False

    def get_sample_data(self) -> Dict:
        """Get sample data from databases for testing"""
        sample_data = {
            "repeat_ids": [],
            "traits": [],
            "chromosomes": [],
            "valid_combinations": [],
        }

        # Get sample repeat IDs from locus database
        if os.path.exists("data/locus_data.db"):
            try:
                conn = sqlite3.connect("data/locus_data.db")
                cursor = conn.cursor()

                # Get sample repeat IDs with their traits
                cursor.execute(
                    """
                    SELECT repeat_id, COALESCE(trait_name, phenotype) as trait, chrom, pos 
                    FROM locus_data 
                    WHERE repeat_id IS NOT NULL 
                    AND (trait_name IS NOT NULL OR phenotype IS NOT NULL)
                    LIMIT 5
                """
                )

                for row in cursor.fetchall():
                    repeat_id, trait, chrom, pos = row
                    sample_data["repeat_ids"].append(repeat_id)
                    sample_data["traits"].append(trait)
                    sample_data["valid_combinations"].append(
                        {
                            "repeat_id": repeat_id,
                            "trait": trait,
                            "chrom": chrom,
                            "pos": pos,
                        }
                    )

                # Remove duplicates from traits
                sample_data["traits"] = list(set(sample_data["traits"]))

                conn.close()
            except Exception as e:
                self.log(
                    f"Warning: Could not get sample data from locus database: {e}",
                    "WARN",
                )

        return sample_data

    def run_basic_tests(self):
        """Run basic endpoint tests"""
        self.log("=== BASIC ENDPOINT TESTS ===")

        # Main pages
        self.test_endpoint("/", description="Home page")
        self.test_endpoint("/about", description="About page")
        self.test_endpoint("/browse_traits", description="Browse traits page")
        self.test_endpoint("/browse_loci", description="Browse loci page")
        self.test_endpoint("/database_status", description="Database status page")
        self.test_endpoint("/setup_instructions", description="Setup instructions page")

        # API endpoints
        self.test_endpoint(
            "/database_status_json", description="Database status JSON API"
        )
        self.test_endpoint("/api/trait_list", description="Trait list API")

    def run_parameterized_tests(self):
        """Run tests with parameters using sample data"""
        self.log("=== PARAMETERIZED TESTS ===")

        sample_data = self.get_sample_data()

        # Test trait overview pages
        if sample_data["traits"]:
            trait = sample_data["traits"][0]
            self.test_endpoint(
                f"/trait_overview/{trait}", description=f"Trait overview for {trait}"
            )

        # Test with valid combinations
        if sample_data["valid_combinations"]:
            combo = sample_data["valid_combinations"][0]

            # Test Manhattan plots with matching trait/repeat_id
            self.test_endpoint(
                "/manhattan_plot",
                params={"trait": combo["trait"], "repeat_id": combo["repeat_id"]},
                description=f"Manhattan plot for {combo['trait']} at {combo['repeat_id']}",
            )

            # Test locus plots
            self.test_endpoint(
                "/locus_plot",
                params={"repeat_id": combo["repeat_id"]},
                description=f"Locus plot for {combo['repeat_id']}",
            )

            # Test with additional parameters
            self.test_endpoint(
                "/locus_plot",
                params={
                    "repeat_id": combo["repeat_id"],
                    "ci_style": "ribbon",
                    "color_by_samples": "true",
                    "count_threshold": "50",
                },
                description=f"Locus plot with custom parameters",
            )

    def run_error_tests(self):
        """Run tests that should return error responses"""
        self.log("=== ERROR HANDLING TESTS ===")

        # Test 404 errors
        self.test_endpoint(
            "/nonexistent_page",
            expected_status=404,
            description="Non-existent page (should return 404)",
        )

        # Test missing parameters
        self.test_endpoint(
            "/manhattan_plot",
            expected_status=400,
            description="Manhattan plot without parameters (should return 400)",
        )
        self.test_endpoint(
            "/locus_plot",
            expected_status=400,
            description="Locus plot without repeat_id (should return 400)",
        )

        # Test invalid parameters
        self.test_endpoint(
            "/manhattan_plot",
            expected_status=404,
            params={"trait": "invalid_trait", "repeat_id": "invalid_id"},
            description="Manhattan plot with invalid parameters",
        )

        self.test_endpoint(
            "/locus_plot",
            expected_status=404,
            params={"repeat_id": "invalid_repeat_id"},
            description="Locus plot with invalid repeat_id",
        )

        # Test invalid trait overview
        self.test_endpoint(
            "/trait_overview/invalid_trait",
            expected_status=404,
            description="Invalid trait overview (should return 404)",
        )

    def run_performance_tests(self):
        """Run performance tests"""
        self.log("=== PERFORMANCE TESTS ===")

        # Test home page multiple times
        times = []
        for i in range(5):
            start = time.time()
            self.test_endpoint("/", description=f"Home page performance test {i+1}")
            times.append((time.time() - start) * 1000)

        avg_time = sum(times) / len(times)
        self.log(f"Average home page response time: {avg_time:.0f}ms")

        if avg_time > 2000:  # 2 seconds
            self.log("Warning: Home page is responding slowly", "WARN")

    def check_server_running(self) -> bool:
        """Check if the Flask server is running"""
        try:
            response = requests.get(f"{self.base_url}/", timeout=5)
            return True
        except requests.exceptions.RequestException:
            return False

    def generate_report(self):
        """Generate a summary report"""
        total_tests = len(self.results)
        successful_tests = sum(1 for r in self.results if r["success"])
        failed_tests = total_tests - successful_tests

        duration = (datetime.now() - self.start_time).total_seconds()

        self.log("=" * 50)
        self.log("TEST SUMMARY REPORT")
        self.log("=" * 50)
        self.log(f"Total tests: {total_tests}")
        self.log(f"Successful: {successful_tests}")
        self.log(f"Failed: {failed_tests}")
        self.log(
            f"Success rate: {(successful_tests/total_tests*100):.1f}%"
            if total_tests > 0
            else "No tests run"
        )
        self.log(f"Total duration: {duration:.1f}s")

        if self.errors:
            self.log("\nERRORS ENCOUNTERED:")
            for i, error in enumerate(self.errors, 1):
                self.log(f"{i}. {error['test']}: {error['error']}")
                if error["response_text"]:
                    self.log(f"   Response preview: {error['response_text'][:200]}...")

        # Performance summary
        if self.results:
            response_times = [
                r["response_time_ms"] for r in self.results if r["success"]
            ]
            if response_times:
                avg_time = sum(response_times) / len(response_times)
                max_time = max(response_times)
                self.log(f"\nPERFORMANCE SUMMARY:")
                self.log(f"Average response time: {avg_time:.0f}ms")
                self.log(f"Slowest response: {max_time:.0f}ms")

        # Save detailed results to file
        self.save_results()

        return failed_tests == 0

    def save_results(self):
        """Save detailed results to a JSON file"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"test_results_{timestamp}.json"

        report_data = {
            "timestamp": self.start_time.isoformat(),
            "base_url": self.base_url,
            "summary": {
                "total_tests": len(self.results),
                "successful_tests": sum(1 for r in self.results if r["success"]),
                "failed_tests": len(self.errors),
            },
            "results": self.results,
            "errors": self.errors,
        }

        try:
            with open(filename, "w") as f:
                json.dump(report_data, f, indent=2)
            self.log(f"Detailed results saved to: {filename}")
        except Exception as e:
            self.log(f"Could not save results file: {e}", "ERROR")


def main():
    """Main test runner"""
    print("STRXplorer Endpoint Testing Script")
    print("=" * 40)

    # Check if server is running
    tester = EndpointTester()

    if not tester.check_server_running():
        print("❌ Flask server is not running!")
        print("Please start the server first:")
        print("   python app.py")
        print("Then run this test script again.")
        sys.exit(1)

    tester.log("Flask server is running, starting tests...")

    # Run all test suites
    try:
        tester.run_basic_tests()
        tester.run_parameterized_tests()
        tester.run_error_tests()
        tester.run_performance_tests()

        # Generate final report
        success = tester.generate_report()

        if success:
            print("\n🎉 All tests passed!")
            sys.exit(0)
        else:
            print("\n⚠️  Some tests failed. Check the errors above.")
            sys.exit(1)

    except KeyboardInterrupt:
        tester.log("Tests interrupted by user", "WARN")
        tester.generate_report()
        sys.exit(1)
    except Exception as e:
        tester.log(f"Unexpected error during testing: {e}", "ERROR")
        sys.exit(1)


if __name__ == "__main__":
    main()
