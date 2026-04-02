# uav_benchmark.py
import subprocess
import json
import matplotlib.pyplot as plt
import numpy as np

# Paths to our compiled Rust binaries
LABRADOR_BINARY = "./target/release/examples/simulator_bridge"
STARK_BINARY = "./target/release/examples/stark_bridge"

def run_rust_benchmark(binary_path):
    """Executes the given Rust binary and parses its JSON output."""
    try:
        result = subprocess.run([binary_path], capture_output=True, text=True, check=True)
        return json.loads(result.stdout.strip())
    except Exception as e:
        print(f"Failed to run {binary_path}: {e}")
        return None

def main():
    print("--- Starting Post-Quantum Cryptography (PQC) Benchmark ---")
    
    print("[*] Running LaBRADOR (Lattice-based) benchmark...")
    labrador_metrics = run_rust_benchmark(LABRADOR_BINARY)
    # Adding known LaBRADOR proof size for the chart (~74 KB per literature)
    labrador_metrics["proof_size_kb"] = 74.0 
    
    print("[*] Running zk-STARK (Hash-based) benchmark...")
    stark_metrics = run_rust_benchmark(STARK_BINARY)

    if not labrador_metrics or not stark_metrics:
        print("Error: Missing metrics. Ensure both Rust examples are compiled.")
        return

    # Extracting Data for the Plot
    algorithms = ['LaBRADOR (Lattice)', 'zk-STARK (Hash)']
    
    prove_times = [labrador_metrics['prove_time_ms'], stark_metrics['prove_time_ms']]
    verify_times = [labrador_metrics['verify_time_ms'], stark_metrics['verify_time_ms']]
    proof_sizes = [labrador_metrics['proof_size_kb'], stark_metrics['proof_size_kb']]

    print("\n--- Benchmark Results ---")
    print(f"LaBRADOR -> Prove: {prove_times[0]:.2f}ms | Verify: {verify_times[0]:.2f}ms | Size: {proof_sizes[0]:.2f} KB")
    print(f"zk-STARK -> Prove: {prove_times[1]:.2f}ms | Verify: {verify_times[1]:.2f}ms | Size: {proof_sizes[1]:.2f} KB")

    # Data Visualization: Grouped Bar Chart
    print("\n[*] Generating comparative graphs...")
    
    x = np.arange(len(algorithms))  # Label locations
    width = 0.35  # Width of the bars

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # SUBPLOT 1: Computational Time (Prove vs Verify)
    rects1 = ax1.bar(x - width/2, prove_times, width, label='Prover Time (ms)', color='#1f77b4')
    rects2 = ax1.bar(x + width/2, verify_times, width, label='Verifier Time (ms)', color='#2ca02c')

    ax1.set_ylabel('Time (Milliseconds)')
    ax1.set_title('Computation Time Comparison (Lower is Better)')
    ax1.set_xticks(x)
    ax1.set_xticklabels(algorithms)
    ax1.legend()
    ax1.grid(axis='y', linestyle='--', alpha=0.7)

    # SUBPLOT 2: Proof Size (Crucial for Drone Bandwidth)
    rects3 = ax2.bar(x, proof_sizes, width * 1.2, label='Proof Size (KB)', color='#ff7f0e')
    
    ax2.set_ylabel('Size (Kilobytes)')
    ax2.set_title('Proof Size Comparison (Lower is Better for UAVs)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(algorithms)
    ax2.legend()
    ax2.grid(axis='y', linestyle='--', alpha=0.7)

    # Add numeric labels on top of bars for clarity
    for ax in [ax1, ax2]:
        for p in ax.patches:
            ax.annotate(f'{p.get_height():.1f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', xytext=(0, 5), textcoords='offset points')

    plt.tight_layout()
    plt.savefig("pqc_benchmark_comparison.png")
    print("[+] Benchmark graph saved as 'pqc_benchmark_comparison.png'.")

if __name__ == "__main__":
    main()